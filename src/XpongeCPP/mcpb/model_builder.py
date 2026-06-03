"""Local model extraction helpers for MCPB workflows."""

from __future__ import annotations

from io import StringIO

from .._core import Molecule, get_template_molecule, has_template, load_mol2 as _core_load_mol2
from .models import MCPBLocalModel
from .selection import _build_atom_to_residue_map, infer_element_symbol


def _single_residue_mol2_text(residue_name, ordered_atoms, bonds):
    lines = [
        "@<TRIPOS>MOLECULE",
        residue_name,
        f"{len(ordered_atoms)} {len(bonds)} 1",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    for atom_index, atom in enumerate(ordered_atoms, start=1):
        lines.append(
            f"{atom_index} {atom['name']} {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f} "
            f"{atom['type']} 1 {residue_name} {atom['charge']:.6f}"
        )
    lines.append("@<TRIPOS>BOND")
    for bond_index, (atom1, atom2) in enumerate(bonds, start=1):
        lines.append(f"{bond_index} {atom1} {atom2} 1")
    return "\n".join(lines) + "\n"


def _residue_intrinsic_bonds(molecule, residue, residue_atom_ids: tuple[int, ...]) -> list[tuple[int, int]]:
    residue_atom_set = set(residue_atom_ids)
    serial_by_atom_id = {atom_id: index + 1 for index, atom_id in enumerate(residue_atom_ids)}
    bonds: list[tuple[int, int]] = []
    for atom1, atom2 in getattr(molecule, "explicit_bonds", None) or []:
        atom1 = int(atom1)
        atom2 = int(atom2)
        if atom1 in residue_atom_set and atom2 in residue_atom_set:
            bonds.append((serial_by_atom_id[atom1], serial_by_atom_id[atom2]))
    if bonds:
        return bonds

    template_name = residue.type_name or residue.name
    if not template_name or not has_template(template_name):
        return bonds
    template = get_template_molecule(template_name)
    if template.residue_count != 1:
        return bonds
    atom_name_to_serial = {
        molecule.atoms[atom_id].name: serial_by_atom_id[atom_id]
        for atom_id in residue_atom_ids
    }
    seen: set[tuple[int, int]] = set()
    for atom1, atom2 in template.explicit_bonds:
        atom1_name = template.atoms[int(atom1)].name
        atom2_name = template.atoms[int(atom2)].name
        serial1 = atom_name_to_serial.get(atom1_name)
        serial2 = atom_name_to_serial.get(atom2_name)
        if serial1 is None or serial2 is None or serial1 == serial2:
            continue
        pair = (serial1, serial2) if serial1 < serial2 else (serial2, serial1)
        if pair in seen:
            continue
        seen.add(pair)
        bonds.append(pair)
    return bonds


def _single_residue_submolecule(molecule, residue_id: int):
    residue = molecule.residues[residue_id]
    residue_atom_ids = tuple(int(atom.index) for atom in residue.atoms)
    ordered_atoms = []
    for atom_id in residue_atom_ids:
        atom = molecule.atoms[atom_id]
        atom_type = atom.type or atom.element or infer_element_symbol(atom.name, residue.name)
        ordered_atoms.append(
            {
                "name": atom.name,
                "type": atom_type,
                "x": atom.x,
                "y": atom.y,
                "z": atom.z,
                "charge": atom.charge,
            }
        )
    bonds = _residue_intrinsic_bonds(molecule, residue, residue_atom_ids)
    text = _single_residue_mol2_text(residue.name, ordered_atoms, bonds)
    return _core_load_mol2(StringIO(text)), residue_atom_ids


def _collect_parent_connect_pairs(molecule) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for atom1, atom2 in getattr(molecule, "residue_links", None) or []:
        atom1 = int(atom1)
        atom2 = int(atom2)
        pair = (atom1, atom2) if atom1 < atom2 else (atom2, atom1)
        if pair in seen:
            continue
        seen.add(pair)
        pairs.append(pair)
    return pairs


def build_local_model(
    molecule,
    residue_ids: tuple[int, ...],
    *,
    name: str,
    extra_connect_pairs: tuple[tuple[int, int], ...] = (),
) -> MCPBLocalModel:
    model = Molecule(name)
    source_atom_ids: list[int] = []
    atom_id_map: dict[int, int] = {}
    residue_id_map: dict[int, int] = {}
    for residue_id in residue_ids:
        residue_model, residue_atom_ids = _single_residue_submolecule(molecule, residue_id)
        atom_offset = model.atom_count
        residue_offset = model.residue_count
        model.add_molecule(residue_model)
        residue_id_map[int(residue_id)] = residue_offset
        for local_index, source_atom_id in enumerate(residue_atom_ids):
            mapped_atom_id = atom_offset + local_index
            atom_id_map[int(source_atom_id)] = mapped_atom_id
            source_atom_ids.append(int(source_atom_id))

    connect_pairs = _collect_parent_connect_pairs(molecule)
    connect_pairs.extend(tuple(pair) for pair in extra_connect_pairs)
    seen_connect: set[tuple[int, int]] = set()
    for atom1, atom2 in connect_pairs:
        if atom1 not in atom_id_map or atom2 not in atom_id_map:
            continue
        mapped = (atom_id_map[atom1], atom_id_map[atom2])
        mapped = mapped if mapped[0] < mapped[1] else (mapped[1], mapped[0])
        if mapped in seen_connect:
            continue
        seen_connect.add(mapped)
        model.add_residue_link(*mapped)

    return MCPBLocalModel(
        molecule=model,
        source_atom_ids=tuple(source_atom_ids),
        atom_id_map=atom_id_map,
        residue_id_map=residue_id_map,
    )


def build_small_and_large_models(request, selection) -> tuple[MCPBLocalModel, MCPBLocalModel]:
    atom_to_residue = _build_atom_to_residue_map(request.molecule)
    ion_residue_ids = {atom_to_residue[atom_id] for atom_id in selection.ion_atom_ids}
    coordinating_residue_ids = {atom_to_residue[atom_id] for atom_id in selection.coordinating_atom_ids}
    small_residue_ids = tuple(sorted(ion_residue_ids | coordinating_residue_ids))
    large_residue_ids = tuple(sorted(set(selection.selected_residue_ids) | set(small_residue_ids)))
    small_model = build_local_model(
        request.molecule,
        small_residue_ids,
        name="MCPB_SMALL",
        extra_connect_pairs=selection.bonded_pairs,
    )
    large_model = build_local_model(
        request.molecule,
        large_residue_ids,
        name="MCPB_LARGE",
        extra_connect_pairs=selection.bonded_pairs,
    )
    return small_model, large_model
