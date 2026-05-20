"""Legacy process-surface wrappers centralized under the compat package."""

from __future__ import annotations

from pathlib import Path

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


def _normalise_link_pair(atom1, atom2):
    if hasattr(atom1, "index") or not isinstance(atom1, (int, str)):
        return _normalise_mol2_bond_pair(atom1, atom2)
    return _normalise_mol2_bond_pair(int(atom1), int(atom2))


def _collect_residue_link_pairs(molecule, residue_links=None):
    connect_pairs = []
    seen = set()
    override = get_legacy_residue_links_override(molecule) or []
    for atom1, atom2 in override:
        pair = _normalise_mol2_bond_pair(atom1, atom2)
        if pair in seen:
            continue
        seen.add(pair)
        connect_pairs.append(pair)
    if residue_links is None:
        residue_links = getattr(molecule, "residue_links", None) or []
    for link in residue_links:
        if hasattr(link, "atom1") and hasattr(link, "atom2"):
            pair = _normalise_link_pair(link.atom1, link.atom2)
        else:
            pair = _normalise_link_pair(link[0], link[1])
        if pair in seen:
            continue
        seen.add(pair)
        connect_pairs.append(pair)
    return connect_pairs


def _patch_saved_pdb_residue_links(molecule, filename, residue_links=None):
    path = Path(filename)
    if not path.exists():
        return None
    with path.open(encoding="utf-8", errors="ignore") as handle:
        lines = handle.read().splitlines()
    end_line = None
    if lines and lines[-1].startswith("END"):
        end_line = lines.pop()
    updated_lines = []
    for line in lines:
        if line.startswith("HETATM"):
            updated_lines.append("ATOM  " + line[6:])
        elif not line.startswith("CONECT"):
            updated_lines.append(line)
    lines = updated_lines

    for atom1, atom2 in _collect_residue_link_pairs(molecule, residue_links=residue_links):
        lines.append(f"CONECT{atom1:5d}{atom2:5d}")
        lines.append(f"CONECT{atom2:5d}{atom1:5d}")
    if end_line is not None:
        lines.append(end_line)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return None


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
    saved_links = None
    if hasattr(target, "residue_links"):
        try:
            saved_links = list(target.residue_links)
        except Exception:
            saved_links = None
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
    if prefix is not None:
        _patch_saved_pdb_residue_links(target, f"{prefix}.pdb", residue_links=saved_links)
    return target


def Save_PDB(molecule, filename, write_cryst1=True):
    target = str(filename)
    save_pdb(molecule, target, write_cryst1)
    _patch_saved_pdb_residue_links(molecule, target)
    return None


def _coerce_link_atom_index(atom):
    if isinstance(atom, (int, str)):
        return int(atom)
    if hasattr(atom, "index") and not callable(atom.index):
        return int(atom.index)
    return int(atom)


def _normalise_mol2_bond_pair(atom1, atom2):
    atom1 = _coerce_link_atom_index(atom1) + 1
    atom2 = _coerce_link_atom_index(atom2) + 1
    if atom1 == atom2:
        raise ValueError("mol2 bond atoms should be different")
    return (atom1, atom2) if atom1 < atom2 else (atom2, atom1)


def _normalise_saved_mol2_bond_pair(atom1, atom2):
    atom1 = int(atom1)
    atom2 = int(atom2)
    if atom1 == atom2:
        raise ValueError("mol2 bond atoms should be different")
    return (atom1, atom2) if atom1 < atom2 else (atom2, atom1)


def _molecule_residue_link_pairs(molecule):
    links = getattr(molecule, "residue_links", None) or []
    pairs = []
    seen = set()
    for link in links:
        if hasattr(link, "atom1") and hasattr(link, "atom2"):
            pair = _normalise_mol2_bond_pair(link.atom1, link.atom2)
        else:
            pair = _normalise_mol2_bond_pair(link[0], link[1])
        if pair in seen:
            continue
        seen.add(pair)
        pairs.append(pair)
    return pairs


def _patch_saved_mol2_residue_links(molecule, filename):
    link_pairs = _molecule_residue_link_pairs(molecule)
    if not link_pairs:
        return None

    with open(filename, encoding="utf-8", errors="ignore") as handle:
        lines = handle.read().splitlines()

    mol2_header = next((i for i, line in enumerate(lines) if line.startswith("@<TRIPOS>MOLECULE")), None)
    bond_header = next((i for i, line in enumerate(lines) if line.startswith("@<TRIPOS>BOND")), None)
    if mol2_header is None or bond_header is None:
        return None

    bond_end = next((i for i in range(bond_header + 1, len(lines)) if lines[i].startswith("@<TRIPOS>")), len(lines))

    existing_pairs = []
    seen = set()
    for line in lines[bond_header + 1:bond_end]:
        words = line.split()
        if len(words) < 4:
            continue
        pair = _normalise_saved_mol2_bond_pair(words[1], words[2])
        if pair in seen:
            continue
        seen.add(pair)
        existing_pairs.append(pair)

    for pair in link_pairs:
        if pair not in seen:
            seen.add(pair)
            existing_pairs.append(pair)

    if len(existing_pairs) == bond_end - bond_header - 1:
        return None

    existing_pairs.sort()
    new_bond_lines = [f"{index:6d} {atom1:5d} {atom2:5d} 1" for index, (atom1, atom2) in enumerate(existing_pairs, start=1)]

    nonempty_after_header = [i for i in range(mol2_header + 1, len(lines)) if lines[i].strip()]
    if len(nonempty_after_header) < 2:
        return None
    count_index = nonempty_after_header[1]
    count_fields = lines[count_index].split()
    if len(count_fields) < 3:
        return None
    status_counts = [int(field) for field in count_fields[3:5]] if len(count_fields) >= 5 else [0, 1]
    lines[count_index] = (
        f"{int(count_fields[0]):6d}{len(existing_pairs):6d}{int(count_fields[2]):6d}"
        f"{status_counts[0]:6d}{status_counts[1]:6d}"
    )

    updated_lines = lines[:bond_header + 1] + new_bond_lines + lines[bond_end:]
    with open(filename, "w", encoding="utf-8") as handle:
        handle.write("\n".join(updated_lines) + "\n")
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
    save_mol2(target, str(filename))
    _patch_saved_mol2_residue_links(target, str(filename))
    return None


def Save_GRO(molecule, filename):
    return save_gro(molecule, str(filename))
