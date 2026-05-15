"""Legacy process-surface wrappers centralized under the compat package."""

from __future__ import annotations

from .._core import (
    Molecule,
    Residue,
    ResidueType,
    add_ions,
    add_molecule,
    add_solvent_box as _core_add_solvent_box,
    get_template_molecule,
    has_template,
    molecule_from_residuetype,
    save_gro,
    save_mol2,
    save_pdb,
    save_sponge_input,
    set_box_padding,
)
from .runtime import get_legacy_residue_links_override


def _single_residue_molecule(value, parameter_name):
    if isinstance(value, Molecule):
        if value.residue_count != 1:
            raise TypeError(f"{parameter_name} molecules should contain exactly one residue")
        return value
    if isinstance(value, ResidueType):
        return molecule_from_residuetype(value)
    if hasattr(value, "name") and has_template(value.name):
        return get_template_molecule(value.name)
    raise TypeError(
        f"{parameter_name} should be a Molecule with one residue, a ResidueType, "
        "or an object whose name matches a registered template"
    )


def Add_Solvent_Box(molecule, solvent, distance, tolerance=2.5, n_solvent=None, seed=0):
    target = molecule if isinstance(molecule, Molecule) else _single_residue_molecule(molecule, "molecule")
    solvent_molecule = solvent if isinstance(solvent, Molecule) else _single_residue_molecule(solvent, "solvent")
    _core_add_solvent_box(target, solvent_molecule, distance, tolerance, n_solvent, seed)
    return target


def Add_Ions(molecule, counts, seed=0, solvent="WAT"):
    normalized = {}
    for key, value in counts.items():
        if isinstance(key, str):
            ion_name = key
        elif isinstance(key, Molecule):
            if key.residue_count != 1:
                raise TypeError("ion molecules should contain exactly one residue")
            ion_name = key.residues[0].name
        elif isinstance(key, ResidueType):
            ion_name = key.name
        elif hasattr(key, "name") and has_template(key.name):
            ion_name = key.name
        else:
            raise TypeError(
                "ion keys should be strings, one-residue Molecules, ResidueType objects, "
                "or template-like objects"
            )
        normalized[str(ion_name)] = int(value)
    return add_ions(molecule, normalized, seed, solvent)


def Add_Molecule(molecule, other):
    return add_molecule(molecule, other)


def Set_Box_Padding(molecule, padding=0.5, center=True):
    return set_box_padding(molecule, padding, center)


def Save_SPONGE_Input(molecule, prefix=None, dirname="."):
    target = molecule
    if isinstance(molecule, Residue):
        target = Molecule(molecule.name)
        target.add_residue(molecule)
    elif isinstance(molecule, ResidueType):
        target = molecule_from_residuetype(molecule)
    elif not isinstance(molecule, Molecule):
        if hasattr(molecule, "name") and has_template(molecule.name):
            target = get_template_molecule(molecule.name)
        else:
            raise TypeError("save_sponge_input expects a Molecule, Residue, ResidueType, or template-like object")
    previous_min_flag = None
    try:
        from ..forcefield.special.min import min_bonded_parameters_enabled

        if min_bonded_parameters_enabled():
            previous_min_flag = False
            target.enable_min_bonded_parameters(True)
    except Exception:
        previous_min_flag = None
    try:
        save_sponge_input(target, "" if prefix is None else str(prefix), str(dirname))
    finally:
        if previous_min_flag is not None:
            target.enable_min_bonded_parameters(False)
    return target


def Save_PDB(molecule, filename, write_cryst1=True):
    target = str(filename)
    save_pdb(molecule, target, write_cryst1)
    override = get_legacy_residue_links_override(molecule)
    if override is None:
        return None
    with open(target, encoding="utf-8", errors="ignore") as handle:
        lines = handle.read().splitlines()
    end_line = None
    if lines and lines[-1].startswith("END"):
        end_line = lines.pop()
    lines = [line for line in lines if not line.startswith("CONECT")]
    for atom1, atom2 in override:
        lines.append(f"CONECT{atom1 + 1:5d}{atom2 + 1:5d}")
        lines.append(f"CONECT{atom2 + 1:5d}{atom1 + 1:5d}")
    if end_line is not None:
        lines.append(end_line)
    with open(target, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return None


def Save_Mol2(molecule, filename=None):
    target = molecule
    if isinstance(molecule, ResidueType):
        target = molecule_from_residuetype(molecule)
    elif not isinstance(molecule, Molecule):
        if hasattr(molecule, "save_as_mol2"):
            return molecule.save_as_mol2(filename)
        if hasattr(molecule, "Save_As_Mol2"):
            return molecule.Save_As_Mol2(filename)
        if hasattr(molecule, "name") and has_template(molecule.name):
            target = get_template_molecule(molecule.name)
        else:
            raise TypeError("Only Molecule, ResidueType, template-like objects, and Assign can be saved as a mol2 file")
    if filename is None:
        name = getattr(target, "name", None) or getattr(molecule, "name", None) or "molecule"
        filename = f"{name}.mol2"
    return save_mol2(target, str(filename))


def Save_GRO(molecule, filename):
    return save_gro(molecule, str(filename))
