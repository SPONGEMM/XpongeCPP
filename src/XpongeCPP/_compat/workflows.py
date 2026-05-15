"""Legacy workflow shims such as gb, min, fep, and build helpers."""

from __future__ import annotations

from collections.abc import Iterable

from .._core import Molecule, Residue, ResidueType
from ..legacy_types import _LegacyResidueTypeHandle

_SAVE_MIN_BONDED_PARAMETERS = False


def set_gb_radius(mol, radius_set="modified_bondi_radii"):
    """Legacy module-level GB helper compatible with Xponge scripts."""
    if hasattr(mol, "set_gb_radius"):
        return mol.set_gb_radius(radius_set)
    if hasattr(mol, "Set_GB_Radius"):
        return mol.Set_GB_Radius(radius_set)
    raise TypeError("set_gb_radius expects a molecule-like object that supports GB radius assignment")


def save_min_bonded_parameters():
    global _SAVE_MIN_BONDED_PARAMETERS
    _SAVE_MIN_BONDED_PARAMETERS = True


def do_not_save_min_bonded_parameters():
    global _SAVE_MIN_BONDED_PARAMETERS
    _SAVE_MIN_BONDED_PARAMETERS = False


def min_bonded_parameters_enabled():
    return _SAVE_MIN_BONDED_PARAMETERS


def build_bonded_force(target, *_args, **_kwargs):
    if isinstance(target, (Molecule, Residue, ResidueType, _LegacyResidueTypeHandle)):
        return target
    raise TypeError(
        "build_bonded_force expects a Molecule, Residue, or ResidueType in "
        "the XpongeCPP first-wave compatibility layer."
    )


def _coerce_legacy_molecule(target):
    from .. import get_template_molecule, has_template, molecule_from_residuetype

    if isinstance(target, Molecule):
        return target
    if isinstance(target, Residue):
        mol = Molecule(target.name)
        mol.add_residue(target)
        return mol
    if isinstance(target, ResidueType):
        return molecule_from_residuetype(target)
    if isinstance(target, _LegacyResidueTypeHandle):
        return get_template_molecule(target.name)
    if hasattr(target, "name") and has_template(target.name):
        return get_template_molecule(target.name)
    raise TypeError(f"The type should be a Molecule, Residue, or ResidueType, but we get {type(target)!r}")


def _coordinate_todo(self, sys_kwarg, ene_kwarg, use_pbc):
    del ene_kwarg
    del use_pbc
    coordinates = self.get_atom_coordinates().tolist()
    if self.box_length is not None:
        box_length = list(self.box_length)
    else:
        coord = self.get_atom_coordinates()
        box_length = list(coord.max(axis=0) - coord.min(axis=0) + 6)
    box_angle = self.box_angle if getattr(self, "box_angle", None) is not None else [90.0, 90.0, 90.0]
    box = [box_length[0], box_length[1], box_length[2], box_angle[0], box_angle[1], box_angle[2]]
    sys_kwarg.setdefault("coordinate", []).append(coordinates)
    sys_kwarg.setdefault("atoms", []).append([atom.name for atom in self.atoms])
    sys_kwarg.setdefault("pbc_box", []).append(box)
    sys_kwarg.setdefault("exclude", []).append([])


_MINDSPONGE_TODOS_INITIALIZED = False


def ensure_mindsponge_todo_support():
    global _MINDSPONGE_TODOS_INITIALIZED
    if _MINDSPONGE_TODOS_INITIALIZED:
        return
    Molecule.Set_MindSponge_Todo("coordinate")(_coordinate_todo)
    # Importing these modules registers first-wave todo handlers.
    from ..forcefield.base import charge_base, mass_base  # noqa: F401

    _MINDSPONGE_TODOS_INITIALIZED = True


def get_mindsponge_system_energy(cls, use_pbc=False):
    ensure_mindsponge_todo_support()
    try:
        from mindsponge import ForceFieldBase
        from mindsponge import Molecule as MindSpongeMolecule
        from mindsponge import set_global_units
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "get_mindsponge_system_energy requires the optional mindsponge/mindspore "
            "dependency stack in the current Python environment."
        ) from exc

    if isinstance(cls, (Molecule, Residue, ResidueType, _LegacyResidueTypeHandle)) or not isinstance(cls, Iterable):
        items = [cls]
    else:
        items = list(cls)

    sys_kwarg = {}
    ene_kwarg = {}
    todos = getattr(Molecule, "_mindsponge_todo", {})
    for item in items:
        mol = _coerce_legacy_molecule(item)
        build_bonded_force(mol)
        for todo in todos.values():
            todo(mol, sys_kwarg, ene_kwarg, use_pbc)

    set_global_units("A", "kcal/mol")
    system = MindSpongeMolecule(**sys_kwarg)
    system.multi_system = len(items)
    energies = []
    for todo in ene_kwarg.values():
        try:
            energies.append(todo["function"](system, ene_kwarg))
        except (TypeError, ValueError) as exc:
            message = str(exc)
            if "NoneType" not in message and "zero dimension" not in message:
                raise

    try:
        energy = ForceFieldBase(energy=energies, exclude_index=sys_kwarg.get("exclude"))
    except ValueError as exc:
        if "zero dimension" not in str(exc):
            raise
        energy = ForceFieldBase(energy=energies)
    return system, energy


Set_GB_Radius = set_gb_radius
Save_Min_Bonded_Parameters = save_min_bonded_parameters
Do_Not_Save_Min_Bonded_Parameters = do_not_save_min_bonded_parameters
