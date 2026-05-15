"""Template registration and residue completion helpers."""

from io import StringIO
from pathlib import Path

import numpy as np

from ._core import (
    Molecule,
    ResidueType,
    get_template_molecule,
    has_template,
    load_mol2 as _core_load_mol2,
    molecule_from_residuetype,
    register_residue_templates_from_mol2_file,
    register_residue_templates_from_mol2_text,
    reorder_atoms_by_template,
    replace_residues,
)
from ._compat.symbols import sync_template_module_globals


def load_mol2(source, ignore_atom_type=False, as_template=False):
    del ignore_atom_type  # First-wave compatibility: current core already tolerates raw MOL2 atom-type strings.
    if isinstance(source, (str, Path)):
        try:
            register_residue_templates_from_mol2_file(str(source))
            sync_template_module_globals()
        except ValueError as exc:
            if as_template or "duplicate atom name in ResidueType" not in str(exc):
                raise
        return _core_load_mol2(str(source))
    if hasattr(source, "read"):
        text = source.read()
        try:
            register_residue_templates_from_mol2_text(text)
            sync_template_module_globals()
        except ValueError as exc:
            if as_template or "duplicate atom name in ResidueType" not in str(exc):
                raise
        return _core_load_mol2(StringIO(text))
    if not as_template:
        return _core_load_mol2(source)
    raise TypeError("load_mol2(...) expects a path-like object or a readable text stream when template registration is required")


def _kabsch(template_positions, fitted_positions):
    template_positions = np.array(template_positions, dtype=np.float32).reshape(-1, 3)
    fitted_positions = np.array(fitted_positions, dtype=np.float32).reshape(-1, 3)
    center1 = np.mean(template_positions, axis=0, keepdims=True)
    center2 = np.mean(fitted_positions, axis=0, keepdims=True)
    if len(template_positions) == 1:
        return np.eye(3), center1.reshape(-1), center2.reshape(-1)
    x_pos = template_positions - center1
    y_pos = fitted_positions - center2
    r_matrix = np.einsum("kj,ki->ij", x_pos, y_pos)
    left, _, right = np.linalg.svd(r_matrix)
    return np.dot(left, right).transpose(), center1.reshape(-1), center2.reshape(-1)


def _template_name_adjacency(template):
    adjacency = {atom.name: [] for atom in template.atoms}
    for atom1, atom2 in template.explicit_bonds:
        atom1_name = template.atoms[int(atom1)].name
        atom2_name = template.atoms[int(atom2)].name
        adjacency[atom1_name].append(atom2_name)
        adjacency[atom2_name].append(atom1_name)
    return adjacency


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


def _aligned_completed_residue(residue):
    template_name = residue.type_name or residue.name
    if not template_name or not has_template(template_name):
        return None
    template = get_template_molecule(template_name)
    if template.residue_count != 1:
        return None
    template_residue = template.residues[0]
    present_names = {atom.name for atom in residue.atoms}
    missing_names = [atom.name for atom in template_residue.atoms if atom.name not in present_names]
    if not missing_names:
        return None
    certified_positions = []
    template_positions = []
    for template_atom in template_residue.atoms:
        if template_atom.name in present_names:
            actual_atom = residue.name2atom(template_atom.name)
            template_positions.append([template_atom.x, template_atom.y, template_atom.z])
            certified_positions.append([actual_atom.x, actual_atom.y, actual_atom.z])
    if not certified_positions:
        return None
    rotation, center1, center2 = _kabsch(template_positions, certified_positions)

    def transform(atom):
        return np.dot(rotation, np.array([atom.x, atom.y, atom.z], dtype=float) - center1) + center2

    adjacency = _template_name_adjacency(template)
    placed_names = set(present_names)
    built_atoms = {
        atom.name: {
            "name": atom.name,
            "type": atom.type,
            "x": atom.x,
            "y": atom.y,
            "z": atom.z,
            "charge": atom.charge,
            "mass": atom.mass,
            "bad_coordinate": atom.bad_coordinate,
            "lj_type_b": atom.lj_type_b,
            "sw_type": atom.sw_type,
            "edip_type": atom.edip_type,
            "gb_radius": atom.gb_radius,
            "gb_scaler": atom.gb_scaler,
            "subsys": atom.subsys,
            "zero_lj_atom": atom.zero_lj_atom,
        }
        for atom in residue.atoms
    }
    unresolved = list(missing_names)
    while unresolved:
        moved = []
        for atom_name in unresolved:
            template_atom = template_residue.name2atom(atom_name)
            transformed_position = transform(template_atom)
            neighbor_names = adjacency.get(atom_name) or [atom.name for atom in template_residue.atoms]
            for neighbor_name in neighbor_names:
                if neighbor_name not in placed_names:
                    continue
                transformed_neighbor = transform(template_residue.name2atom(neighbor_name))
                anchored_neighbor = built_atoms[neighbor_name]
                anchored = transformed_position - transformed_neighbor + np.array(
                    [anchored_neighbor["x"], anchored_neighbor["y"], anchored_neighbor["z"]],
                    dtype=float,
                )
                built_atoms[atom_name] = {
                    "name": template_atom.name,
                    "type": template_atom.type,
                    "x": float(anchored[0]),
                    "y": float(anchored[1]),
                    "z": float(anchored[2]),
                    "charge": template_atom.charge,
                    "mass": template_atom.mass,
                    "bad_coordinate": True,
                    "lj_type_b": template_atom.lj_type_b,
                    "sw_type": template_atom.sw_type,
                    "edip_type": template_atom.edip_type,
                    "gb_radius": template_atom.gb_radius,
                    "gb_scaler": template_atom.gb_scaler,
                    "subsys": template_atom.subsys,
                    "zero_lj_atom": template_atom.zero_lj_atom,
                }
                placed_names.add(atom_name)
                moved.append(atom_name)
                break
        if not moved:
            for atom_name in unresolved:
                anchored = transform(template_residue.name2atom(atom_name))
                template_atom = template_residue.name2atom(atom_name)
                built_atoms[atom_name] = {
                    "name": template_atom.name,
                    "type": template_atom.type,
                    "x": float(anchored[0]),
                    "y": float(anchored[1]),
                    "z": float(anchored[2]),
                    "charge": template_atom.charge,
                    "mass": template_atom.mass,
                    "bad_coordinate": True,
                    "lj_type_b": template_atom.lj_type_b,
                    "sw_type": template_atom.sw_type,
                    "edip_type": template_atom.edip_type,
                    "gb_radius": template_atom.gb_radius,
                    "gb_scaler": template_atom.gb_scaler,
                    "subsys": template_atom.subsys,
                    "zero_lj_atom": template_atom.zero_lj_atom,
                }
            break
        unresolved = [name for name in unresolved if name not in moved]

    ordered_names = [atom.name for atom in residue.atoms]
    ordered_names.extend(name for name in missing_names if name not in present_names)
    ordered_atoms = [built_atoms[name] for name in ordered_names]
    serial_by_name = {name: index + 1 for index, name in enumerate(ordered_names)}
    bonds = []
    for atom1, atom2 in template.explicit_bonds:
        atom1_name = template.atoms[int(atom1)].name
        atom2_name = template.atoms[int(atom2)].name
        if atom1_name in serial_by_name and atom2_name in serial_by_name:
            bonds.append((serial_by_name[atom1_name], serial_by_name[atom2_name]))
    replacement = _core_load_mol2(StringIO(_single_residue_mol2_text(template_name, ordered_atoms, bonds)))
    replacement_residue = replacement.residues[0]
    for atom in replacement_residue.atoms:
        atom_data = built_atoms[atom.name]
        atom.type = atom_data["type"]
        atom.x = atom_data["x"]
        atom.y = atom_data["y"]
        atom.z = atom_data["z"]
        atom.charge = atom_data["charge"]
        atom.mass = atom_data["mass"]
        atom.bad_coordinate = atom_data["bad_coordinate"]
        atom.lj_type_b = atom_data["lj_type_b"]
        atom.sw_type = atom_data["sw_type"]
        atom.edip_type = atom_data["edip_type"]
        atom.gb_radius = atom_data["gb_radius"]
        atom.gb_scaler = atom_data["gb_scaler"]
        atom.subsys = atom_data["subsys"]
        atom.zero_lj_atom = atom_data["zero_lj_atom"]
    return replacement


def _legacy_add_missing_atoms(self):
    replacements = {}
    for residue in self.residues:
        replacement = _aligned_completed_residue(residue)
        if replacement is not None:
            replacements[int(residue.index)] = replacement
    if replacements:
        replaced_residue_ids = set(replacements)
        atom_map_before = {
            int(atom.index): (int(residue.index), atom.name)
            for residue in self.residues
            for atom in residue.atoms
        }
        preserved_links = []
        for atom1, atom2 in self.residue_links:
            descriptor1 = atom_map_before[int(atom1)]
            descriptor2 = atom_map_before[int(atom2)]
            if descriptor1[0] in replaced_residue_ids or descriptor2[0] in replaced_residue_ids:
                preserved_links.append(tuple(sorted((descriptor1, descriptor2))))
        replace_residues(self, replacements, [], False)
        atom_map_after = {
            int(atom.index): (int(residue.index), atom.name)
            for residue in self.residues
            for atom in residue.atoms
        }
        current_links = {
            tuple(sorted((atom_map_after[int(atom1)], atom_map_after[int(atom2)])))
            for atom1, atom2 in self.residue_links
        }
        for descriptor1, descriptor2 in preserved_links:
            normalized = tuple(sorted((descriptor1, descriptor2)))
            if normalized in current_links:
                continue
            residue1 = self.residues[descriptor1[0]]
            residue2 = self.residues[descriptor2[0]]
            self.add_residue_link(residue1.name2atom(descriptor1[1]), residue2.name2atom(descriptor2[1]))
            current_links.add(normalized)
    return self
