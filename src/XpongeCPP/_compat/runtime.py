"""Centralized installation of legacy runtime patches."""

from __future__ import annotations

from .imports import Xdict
from ..legacy_types import _legacy_residue_links_override
from .._core import Assign, Molecule, Residue, ResidueType
from ..legacy_types import (
    _AtomIndexProxy,
    _legacy_add_residue,
    _legacy_add_residue_link,
    _legacy_add_residue_links,
    _legacy_clear_residue_links,
    _legacy_get_residuetype,
    _legacy_get_residue_links,
    _legacy_get_residue_links_copy,
    _legacy_make_residue_like,
    _legacy_set_residue_links,
)
from .._core import molecule_from_residuetype
from ..template_ops import _legacy_add_missing_atoms


def _legacy_residue_factory(template_like, directly_copy=True):
    return _legacy_make_residue_like(template_like, directly_copy=directly_copy)


_core_residuetype_add_atom = ResidueType.add_atom


def _legacy_residuetype_add_atom(self, name, atom_type, x, y, z, charge=0.0, mass=0.0):
    if hasattr(atom_type, "name"):
        atom_type = atom_type.name
    return _core_residuetype_add_atom(self, name, str(atom_type), x, y, z, charge, mass)


def _legacy_residuetype_atoms(self):
    return molecule_from_residuetype(self).residues[0].atoms


class _LegacyCallableTypeProxy:
    """Thin callable proxy that accepts legacy keyword constructors.

    Old Xponge scripts often call core classes as ``ResidueType(name="ALA")`` or
    ``Molecule(name="foo")``.  The pybind-backed XpongeCPP classes mostly expect
    positional string arguments.  This proxy keeps the class attribute surface
    while translating the high-frequency keyword form.
    """

    def __init__(self, target, keyword_name="name"):
        self._target = target
        self._keyword_name = keyword_name

    def __call__(self, *args, **kwargs):
        if kwargs:
            if not args and set(kwargs) == {self._keyword_name}:
                return self._target(kwargs[self._keyword_name])
            return self._target(*args, **kwargs)
        return self._target(*args)

    def __getattr__(self, item):
        return getattr(self._target, item)

    def __repr__(self):
        return repr(self._target)


def get_legacy_residue_links_override(molecule):
    """Return the Python-side residue-link override list when present."""

    return _legacy_residue_links_override.get(molecule)


def _legacy_set_save_sponge_input(cls, keyname):
    def wrapper(func):
        cls._save_functions[str(keyname)] = func
        return func

    return wrapper


def _legacy_del_save_sponge_input(cls, keyname):
    cls._save_functions.pop(str(keyname), None)


def _legacy_set_mindsponge_todo(cls, keyname):
    def wrapper(func):
        cls._mindsponge_todo[str(keyname)] = func
        return func

    return wrapper


def _legacy_del_mindsponge_todo(cls, keyname):
    cls._mindsponge_todo.pop(str(keyname), None)


def install_legacy_runtime_patches(namespace: dict | None = None):
    """Install legacy object-surface shims on core runtime classes.

    This keeps class monkeypatches and module-level ``Residue(...)`` factory
    wiring in one place instead of scattering them through ``__init__.py``.
    """

    if not hasattr(Molecule, "_save_functions"):
        Molecule._save_functions = Xdict()
    if not hasattr(Molecule, "_mindsponge_todo"):
        Molecule._mindsponge_todo = Xdict()

    ResidueType.get_type = staticmethod(_legacy_get_residuetype)
    ResidueType.Get_Type = staticmethod(_legacy_get_residuetype)
    ResidueType.add_atom = _legacy_residuetype_add_atom
    ResidueType.addAtom = _legacy_residuetype_add_atom
    ResidueType.Add_Atom = _legacy_residuetype_add_atom
    ResidueType.atoms = property(_legacy_residuetype_atoms)

    Molecule.add_residue = _legacy_add_residue
    Molecule.Add_Residue = _legacy_add_residue
    Molecule.add_residue_link = _legacy_add_residue_link
    Molecule.Add_Residue_Link = _legacy_add_residue_link
    Molecule.add_residue_links = _legacy_add_residue_links
    Molecule.Add_Residue_Links = _legacy_add_residue_links
    Molecule.clear_residue_links = _legacy_clear_residue_links
    Molecule.Clear_Residue_Links = _legacy_clear_residue_links
    Molecule.set_residue_links = _legacy_set_residue_links
    Molecule.Set_Residue_Links = _legacy_set_residue_links
    Molecule.get_residue_links = _legacy_get_residue_links_copy
    Molecule.Get_Residue_Links = _legacy_get_residue_links_copy
    Molecule.residue_links = property(_legacy_get_residue_links, _legacy_set_residue_links)
    Molecule.atom_index = property(lambda self: _AtomIndexProxy())
    Molecule.add_missing_atoms = _legacy_add_missing_atoms
    Molecule.Add_Missing_Atoms = _legacy_add_missing_atoms
    Molecule.set_save_sponge_input = classmethod(_legacy_set_save_sponge_input)
    Molecule.Set_Save_SPONGE_Input = Molecule.set_save_sponge_input
    Molecule.del_save_sponge_input = classmethod(_legacy_del_save_sponge_input)
    Molecule.Del_Save_SPONGE_Input = Molecule.del_save_sponge_input
    Molecule.set_mindsponge_todo = classmethod(_legacy_set_mindsponge_todo)
    Molecule.Set_MindSponge_Todo = Molecule.set_mindsponge_todo
    Molecule.del_mindsponge_todo = classmethod(_legacy_del_mindsponge_todo)
    Molecule.Del_MindSponge_Todo = Molecule.del_mindsponge_todo

    Residue.unterminal = lambda self: self
    Residue.Unterminal = Residue.unterminal
    Residue.UnTerminal = Residue.unterminal

    Assign.addAtom = Assign.add_atom
    Assign.Add_Bond = Assign.add_bond
    Assign.addBond = Assign.add_bond
    Assign.deleteAtom = Assign.delete_atom
    Assign.Delete_Atom = Assign.delete_atom
    Assign.deleteBond = Assign.delete_bond
    Assign.Delete_Bond = Assign.delete_bond

    if namespace is not None:
        namespace["Residue"] = _legacy_residue_factory
        namespace["ResidueType"] = _LegacyCallableTypeProxy(ResidueType)
        namespace["Molecule"] = _LegacyCallableTypeProxy(Molecule)

    return _legacy_residue_factory
