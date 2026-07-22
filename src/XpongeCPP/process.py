"""Geometry and molecule processing helpers built on top of the core facade."""

import math
import os
import shutil
import subprocess
import tempfile
from importlib import import_module

import numpy as np

from ._core import (
    Molecule,
    Residue,
    ResidueType,
    get_template_molecule,
    has_template,
    load_coordinate,
    reorder_atoms_by_template,
    replace_residues,
)
from ._compat.process import (
    Add_Ions,
    Add_Molecule,
    Add_Solvent_Box,
    Save_GRO,
    Save_Mol2,
    Save_PDB,
    Save_SPONGE_Input,
    Set_Box_Padding,
    _single_residue_molecule,
    save_sponge_input_raw,
)


add_solvent_box = Add_Solvent_Box
Save_Sponge_Input = Save_SPONGE_Input


def _rotation_matrix(axis, angle):
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm == 0.0:
        return np.eye(3)
    x, y, z = axis / norm
    cosine = math.cos(angle)
    sine = math.sin(angle)
    one_minus_cosine = 1.0 - cosine
    return np.array([
        [cosine + x * x * one_minus_cosine, x * y * one_minus_cosine - z * sine, x * z * one_minus_cosine + y * sine],
        [y * x * one_minus_cosine + z * sine, cosine + y * y * one_minus_cosine, y * z * one_minus_cosine - x * sine],
        [z * x * one_minus_cosine - y * sine, z * y * one_minus_cosine + x * sine, cosine + z * z * one_minus_cosine],
    ]).transpose()


def _molecule_coordinates(molecule):
    return np.array([[atom.x, atom.y, atom.z] for atom in molecule.atoms], dtype=float)


def _set_molecule_coordinates(molecule, coordinates):
    for atom, coordinate in zip(molecule.atoms, coordinates):
        atom.x = float(coordinate[0])
        atom.y = float(coordinate[1])
        atom.z = float(coordinate[2])


def _bond_adjacency(molecule):
    adjacency = {int(atom.index): set() for atom in molecule.atoms}
    for atom1, atom2 in molecule.explicit_bonds:
        atom1 = int(atom1)
        atom2 = int(atom2)
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)
    for atom1, atom2 in molecule.residue_links:
        atom1 = int(atom1)
        atom2 = int(atom2)
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)
    return adjacency


def _divide_into_two_parts(molecule, atom1, atom2):
    atom1_index = int(atom1.index)
    atom2_index = int(atom2.index)
    adjacency = _bond_adjacency(molecule)
    adjacency[atom1_index].discard(atom2_index)
    adjacency[atom2_index].discard(atom1_index)
    left = set()
    stack = [atom1_index]
    while stack:
        current = stack.pop()
        if current in left:
            continue
        left.add(current)
        stack.extend(adjacency[current] - left)
    if atom2_index in left:
        raise ValueError("atom1 and atom2 remain connected after splitting the bond")
    right = set()
    stack = [atom2_index]
    while stack:
        current = stack.pop()
        if current in right:
            continue
        right.add(current)
        stack.extend(adjacency[current] - right)
    return sorted(left), sorted(right)


def impose_bond(molecule, atom1, atom2, length):
    coordinates = _molecule_coordinates(molecule)
    _, atom2_friends = _divide_into_two_parts(molecule, atom1, atom2)
    r0 = coordinates[int(atom2.index)] - coordinates[int(atom1.index)]
    l0 = np.linalg.norm(r0)
    if l0 == 0.0:
        coordinates[int(atom2.index)] += math.sqrt(1.0 / 3.0)
        r0 = coordinates[int(atom2.index)] - coordinates[int(atom1.index)]
        l0 = np.linalg.norm(r0)
    dr = (float(length) / l0 - 1.0) * r0
    coordinates[atom2_friends] += dr
    _set_molecule_coordinates(molecule, coordinates)


def impose_angle(molecule, atom1, atom2, atom3, angle):
    coordinates = _molecule_coordinates(molecule)
    _, atom3_friends = _divide_into_two_parts(molecule, atom2, atom3)
    r12 = coordinates[int(atom1.index)] - coordinates[int(atom2.index)]
    r23 = coordinates[int(atom3.index)] - coordinates[int(atom2.index)]
    angle0 = math.acos(float(np.dot(r12, r23) / np.linalg.norm(r23) / np.linalg.norm(r12)))
    delta_angle = float(angle) - angle0
    rotation = _rotation_matrix(np.cross(r12, r23), delta_angle)
    origin = coordinates[int(atom2.index)]
    coordinates[atom3_friends] = np.dot(coordinates[atom3_friends] - origin, rotation) + origin
    _set_molecule_coordinates(molecule, coordinates)


def impose_dihedral(molecule, atom1, atom2, atom3, atom4, dihedral, keep_atom3=False):
    coordinates = _molecule_coordinates(molecule)
    if not keep_atom3:
        _, rotate_friends = _divide_into_two_parts(molecule, atom2, atom3)
    else:
        _, rotate_friends = _divide_into_two_parts(molecule, atom3, atom4)
    r12 = coordinates[int(atom1.index)] - coordinates[int(atom2.index)]
    r23 = coordinates[int(atom3.index)] - coordinates[int(atom2.index)]
    r34 = coordinates[int(atom3.index)] - coordinates[int(atom4.index)]
    r12xr23 = np.cross(r12, r23)
    r34xr23 = np.cross(r34, r23)
    cosine = float(np.dot(r12xr23, r34xr23) / np.linalg.norm(r12xr23) / np.linalg.norm(r34xr23))
    cosine = max(-0.999999, min(cosine, 0.999999))
    dihedral0 = math.acos(cosine)
    dihedral0 = math.pi - math.copysign(dihedral0, float(np.dot(np.cross(r34xr23, r12xr23), r23)))
    delta_angle = float(dihedral) - dihedral0
    rotation = _rotation_matrix(r23, delta_angle)
    origin = coordinates[int(atom3.index)]
    coordinates[rotate_friends] = np.dot(coordinates[rotate_friends] - origin, rotation) + origin
    _set_molecule_coordinates(molecule, coordinates)


def main_axis_rotate(molecule, direction_long=None, direction_middle=None, direction_short=None):
    direction_long = np.array(direction_long if direction_long is not None else [0, 0, 1], dtype=float)
    direction_middle = np.array(direction_middle if direction_middle is not None else [0, 1, 0], dtype=float)
    direction_short = np.array(direction_short if direction_short is not None else [1, 0, 0], dtype=float)
    coordinates = _molecule_coordinates(molecule)
    center = np.zeros(3, dtype=float)
    total_mass = 0.0
    for atom, coordinate in zip(molecule.atoms, coordinates):
        total_mass += atom.mass
        center += atom.mass * coordinate
    center /= total_mass
    inertia = np.zeros((3, 3), dtype=float)
    for atom, coordinate in zip(molecule.atoms, coordinates):
        x, y, z = coordinate - center
        inertia += atom.mass * np.array([
            [y * y + z * z, -x * y, -x * z],
            [-x * y, x * x + z * z, -y * z],
            [-x * z, -y * z, x * x + y * y],
        ])
    eigenvalues, eigenvectors = np.linalg.eig(inertia)
    order = np.argsort(eigenvalues)
    target = np.vstack([direction_short, direction_middle, direction_long])
    source = np.vstack((eigenvectors[:, order[2]], eigenvectors[:, order[1]], eigenvectors[:, order[0]]))
    rotation = np.dot(target, np.linalg.inv(source))
    rotated = np.dot(coordinates - center, rotation) + center
    _set_molecule_coordinates(molecule, rotated)


def get_peptide_from_sequence(sequence, charged_terminal=True):
    if not isinstance(sequence, str) or len(sequence) <= 1:
        raise AssertionError("sequence should be a string longer than 1")
    mapping = {
        "A": "ALA", "G": "GLY", "V": "VAL", "L": "LEU", "I": "ILE", "P": "PRO",
        "F": "PHE", "Y": "TYR", "W": "TRP", "S": "SER", "T": "THR", "C": "CYS",
        "M": "MET", "N": "ASN", "Q": "GLN", "D": "ASP", "E": "GLU", "K": "LYS",
        "R": "ARG", "H": "HIS",
    }
    head = mapping[sequence[0]]
    tail = mapping[sequence[-1]]
    if charged_terminal:
        head = "N" + head
        tail = "C" + tail
    result = get_template_molecule(head)
    for code in sequence[1:-1]:
        result = result + get_template_molecule(mapping[code])
    result = result + get_template_molecule(tail)
    return result


def h_mass_repartition(molecule, repartition_mass=1.1, repartition_rate=3, exclude_residue_name="WAT"):
    adjacency = _bond_adjacency(molecule)
    for residue in molecule.residues:
        if residue.name == exclude_residue_name:
            continue
        template_name = residue.type_name or residue.name
        template_neighbors = None
        if has_template(template_name):
            template = get_template_molecule(template_name)
            template_neighbors = {atom.name: set() for atom in template.residues[0].atoms}
            for atom1, atom2 in template.explicit_bonds:
                atom1_name = template.atoms[int(atom1)].name
                atom2_name = template.atoms[int(atom2)].name
                template_neighbors[atom1_name].add(atom2_name)
                template_neighbors[atom2_name].add(atom1_name)
        residue_indices = {int(atom.index) for atom in residue.atoms}
        for atom in residue.atoms:
            if atom.mass > repartition_mass:
                continue
            heavy_atom = None
            if template_neighbors is not None:
                neighbor_names = sorted(template_neighbors.get(atom.name, set()))
                if len(neighbor_names) == 1:
                    heavy_atom = residue.name2atom(neighbor_names[0])
            if heavy_atom is None:
                neighbors = [index for index in adjacency[int(atom.index)] if index in residue_indices]
                if len(neighbors) != 1:
                    continue
                heavy_atom = molecule.atoms[neighbors[0]]
            origin_mass = atom.mass
            atom.mass *= repartition_rate
            heavy_atom.mass -= atom.mass - origin_mass


def solvent_replace(molecule, select, toreplace, sort=True):
    if not callable(select):
        if isinstance(select, Molecule):
            if select.residue_count != 1:
                raise TypeError("select molecules should contain exactly one residue")
            resname = select.residues[0].name
        elif isinstance(select, ResidueType):
            resname = select.name
        elif hasattr(select, "name"):
            resname = select.name
        else:
            raise TypeError("select should be callable, a Molecule, a ResidueType, or an object with a name")
        select = lambda residue, target_name=resname: residue.name == target_name

    selected = []
    residue_sort_keys = []
    for residue in molecule.residues:
        if select(residue):
            selected.append(int(residue.index))
            residue_sort_keys.append(float("inf"))
        else:
            residue_sort_keys.append(float("-inf"))

    np.random.shuffle(selected)
    replacements = {}
    count = 0
    for key, value in toreplace.items():
        replacement = _single_residue_molecule(key, "toreplace keys")
        value = int(value)
        indices = selected[count:count + value]
        count += value
        for residue_index in indices:
            replacements[residue_index] = replacement
            residue_sort_keys[residue_index] = float(count)
    replace_residues(molecule, replacements, residue_sort_keys, sort)


def sort_atoms_by(mol, template):
    if not isinstance(mol, Molecule) or not isinstance(template, Molecule):
        raise TypeError("the type of the input should be Xponge.Molecule")
    reorder_atoms_by_template(mol, template)


def optimize(mol, step=2000, only_bad_coordinate=True, dt=1e-8, pbc=True, extra_commands=None):
    executable = "SPONGE" if pbc else "SPONGE_NOPBC"
    if shutil.which(executable) is None:
        raise RuntimeError(f"{executable} executable is required for optimize()")
    if extra_commands is not None and not hasattr(extra_commands, "items"):
        raise TypeError("extra_commands should be a mapping of command names to values")
    with tempfile.TemporaryDirectory() as tempdir:
        prefix = "temp"
        save_prefix = os.path.join(tempdir, prefix)
        output_prefix = os.path.join(tempdir, "min")
        mol.enable_min_bonded_parameters(True)
        box_backup = None
        if not pbc:
            box_backup = (list(mol.box_length), list(mol.box_angle))
            mol.box_length = [999.0, 999.0, 999.0]
            mol.box_angle = [90.0, 90.0, 90.0]
        try:
            import_module("XpongeCPP").Save_SPONGE_Input(mol, prefix=prefix, dirname=tempdir)
        finally:
            mol.enable_min_bonded_parameters(False)
            if box_backup is not None:
                mol.box_length = box_backup[0]
                mol.box_angle = box_backup[1]
        mdin_name = os.path.join(tempdir, "mdin.txt")
        with open(mdin_name, "w", encoding="utf-8") as handle:
            handle.write(
                f"""temp
default_in_file_prefix = {save_prefix}
rst = {output_prefix}
crd = {save_prefix}.dat
box = {save_prefix}.box
mdout = {output_prefix}.out
mdinfo = {output_prefix}.info
mode = minimization
minimization_dynamic_dt = 1
step_limit = {step}
write_information_interval = {step}
molecule_map_output = 1
dont_check_input = 1
"""
            )
            if extra_commands:
                for command, value in extra_commands.items():
                    handle.write(f"{command} = {value}\n")
        command = [executable, "-mdin", mdin_name, "-dt", str(dt)]
        if only_bad_coordinate:
            command.extend(["-mass_in_file", f"{save_prefix}_fake_mass.txt"])
        result = subprocess.run(command, check=False, text=True, capture_output=True)
        if result.returncode != 0:
            message = result.stderr.strip() or result.stdout.strip() or "The optimization failed"
            raise RuntimeError(message)
        coordinate_filename = f"{output_prefix}_coordinate.txt"
        if not os.path.exists(coordinate_filename):
            raise RuntimeError("The optimization did not produce coordinate output")
        load_coordinate(coordinate_filename, mol)


class Region:
    def __init__(self, side="in", boundary=False):
        self._side = True
        self.side = side
        self.boundary = boundary

    @property
    def side(self):
        return "in" if self._side else "out"

    @side.setter
    def side(self, side):
        if side == "in":
            self._side = True
        elif side == "out":
            self._side = False
        else:
            raise ValueError("side should be 'in' or 'out'")


class IntersectRegion:
    def __init__(self, *regions):
        self.regions = regions

    def __contains__(self, item):
        return all(item in region for region in self.regions)


class UnionRegion:
    def __init__(self, *regions):
        self.regions = regions

    def __contains__(self, item):
        return any(item in region for region in self.regions)


class BlockRegion(Region):
    def __init__(self, x_low, y_low, z_low, x_high, y_high, z_high, side="in", boundary=False):
        self.x_low = x_low
        self.y_low = y_low
        self.z_low = z_low
        self.x_high = x_high
        self.y_high = y_high
        self.z_high = z_high
        super().__init__(side, boundary)

    def __contains__(self, item):
        if self.boundary:
            ans = self.x_low <= item[0] <= self.x_high and self.y_low <= item[1] <= self.y_high and self.z_low <= item[2] <= self.z_high
        else:
            ans = self.x_low < item[0] < self.x_high and self.y_low < item[1] < self.y_high and self.z_low < item[2] < self.z_high
        return ans if self._side else not ans


class SphereRegion(Region):
    def __init__(self, x, y, z, r, side="in", boundary=False):
        self.x = x
        self.y = y
        self.z = z
        self._r2 = r * r
        super().__init__(side, boundary)

    def __contains__(self, item):
        ans = (item[0] - self.x) ** 2 + (item[1] - self.y) ** 2 + (item[2] - self.z) ** 2
        ans = ans <= self._r2 if self.boundary else ans < self._r2
        return ans if self._side else not ans


class FrustumRegion(Region):
    def __init__(self, x1, y1, z1, r1, x2, y2, z2, r2, side="in", boundary=False):
        self.r1 = r1
        self.r2 = r2
        self.o1 = np.array([x1, y1, z1], dtype=np.float32)
        self.o2 = np.array([x2, y2, z2], dtype=np.float32)
        self.axis = self.o2 - self.o1
        self.length = np.linalg.norm(self.axis)
        self.axis /= self.length
        self.k = (r2 - r1) / self.length
        super().__init__(side, boundary)

    def __contains__(self, item):
        coordinate = np.array(item) - self.o1
        length = np.linalg.norm(coordinate)
        projection = np.dot(coordinate, self.axis)
        distance = length * length - projection * projection
        radius = self.r1 + self.k * projection
        if self.boundary:
            ans = self.length >= projection >= 0 and distance <= radius * radius
        else:
            ans = self.length > projection > 0 and distance < radius * radius
        return ans if self._side else not ans


class PrismRegion(Region):
    def __init__(self, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, side="in", boundary=False):
        self.l0 = np.array([x0, y0, z0], dtype=np.float32)
        self.l1 = np.array([x1, y1, z1], dtype=np.float32)
        self.l2 = np.array([x2, y2, z2], dtype=np.float32)
        self.l3 = np.array([x3, y3, z3], dtype=np.float32)
        self.n3 = np.cross(self.l1, self.l2)
        self.n3 /= np.linalg.norm(self.n3)
        self.n2 = np.cross(self.l3, self.l1)
        self.n2 /= np.linalg.norm(self.n2)
        self.n1 = np.cross(self.l2, self.l3)
        self.n1 /= np.linalg.norm(self.n1)
        self.length = np.array([np.dot(self.l1, self.n1), np.dot(self.l2, self.n2), np.dot(self.l3, self.n3)])
        assert np.all(self.length > 0), "The basis vectors should mmet the right-handed axis system requirements"
        super().__init__(side, boundary)

    def __contains__(self, item):
        coordinate = np.array(item) - self.l0
        if self.boundary:
            ans = 0 <= np.dot(coordinate, self.n1) <= self.length[0] and 0 <= np.dot(coordinate, self.n2) <= self.length[1] and 0 <= np.dot(coordinate, self.n3) <= self.length[2]
        else:
            ans = 0 < np.dot(coordinate, self.n1) < self.length[0] and 0 < np.dot(coordinate, self.n2) < self.length[1] and 0 < np.dot(coordinate, self.n3) < self.length[2]
        return ans if self._side else not ans


def _length_angle_from_basis_vectors(l1, l2, l3):
    l1 = np.array(l1, dtype=float)
    l2 = np.array(l2, dtype=float)
    l3 = np.array(l3, dtype=float)
    la = np.linalg.norm(l1)
    lb = np.linalg.norm(l2)
    lc = np.linalg.norm(l3)
    alpha = math.degrees(math.acos(np.dot(l2, l3) / lb / lc))
    beta = math.degrees(math.acos(np.dot(l1, l3) / la / lc))
    gamma = math.degrees(math.acos(np.dot(l1, l2) / la / lb))
    return [la, lb, lc], [alpha, beta, gamma]


def _basis_vectors_from_length_and_angle(la, lb, lc, alpha, beta, gamma):
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)
    ax, ay, az = la, 0.0, 0.0
    bx, by, bz = lb * math.cos(gamma), lb * math.sin(gamma), 0.0
    cx = lc * math.cos(beta)
    cy = lc * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz2 = lc * lc - cx * cx - cy * cy
    cz = math.sqrt(max(cz2, 0.0))
    return [(ax, ay, az), (bx, by, bz), (cx, cy, cz)]


class Lattice:
    styles = {}

    def __init__(
        self,
        style="custom",
        basis_molecule=None,
        scale=None,
        origin=None,
        cell_length=None,
        cell_angle=None,
        basis_position=None,
        spacing=None,
        periodic_bonds=None,
        periodic_cutoff=3,
    ):
        self.basis_molecule = basis_molecule
        if style == "custom" or style.startswith("template:"):
            self.scale = 1 if scale is None else scale
            self.origin = [0, 0, 0] if origin is None else origin
            self.cell_length = [1, 1, 1] if cell_length is None else cell_length
            self.cell_angle = np.array([90, 90, 90] if cell_angle is None else cell_angle, dtype=float)
            assert np.all((self.cell_angle > 0) & (self.cell_angle < 180)), "the cell angle should be in the range (0, 180)"
            self.spacing = [0, 0, 0] if spacing is None else spacing
            self.basis_position = [] if basis_position is None else basis_position
        else:
            old_style = Lattice.styles[style]
            self.scale = scale
            self.origin = old_style.origin
            self.cell_length = old_style.cell_length
            self.cell_angle = old_style.cell_angle
            self.basis_position = old_style.basis_position
            self.spacing = old_style.spacing
        if not style.startswith("template:") and self.basis_molecule is None:
            raise ValueError("basis molecule should not be None for a non-template lattice")
        if not style.startswith("template:") and self.scale is None:
            raise ValueError("scale should not be None for a non-template lattice")
        if style.startswith("template:"):
            Lattice.styles[style.split(":")[1].strip()] = self
        self.periodic_bonds = periodic_bonds
        self.current_unbonded_periodic_atoms = set()
        self.periodic_cutoff = periodic_cutoff * periodic_cutoff

    def _process_periodic_bonds(self, mol, residue, box):
        if not self.periodic_bonds:
            return
        for name1, name2 in self.periodic_bonds:
            self.current_unbonded_periodic_atoms.add((residue.name2atom(name1), name2))
            self.current_unbonded_periodic_atoms.add((residue.name2atom(name2), name1))
        remove = set()
        for atom, name in self.current_unbonded_periodic_atoms:
            atom2 = residue.name2atom(name)
            dx = atom.x - atom2.x
            dx -= math.floor(dx / (box.x_high - box.x_low) + 0.5) * (box.x_high - box.x_low)
            dy = atom.y - atom2.y
            dy -= math.floor(dy / (box.y_high - box.y_low) + 0.5) * (box.y_high - box.y_low)
            dz = atom.z - atom2.z
            dz -= math.floor(dz / (box.z_high - box.z_low) + 0.5) * (box.z_high - box.z_low)
            if dx * dx + dy * dy + dz * dz < self.periodic_cutoff:
                mol.add_residue_link(int(atom.index), int(atom2.index))
                remove.add((atom, name))
                remove.add((atom2, atom.name))
        self.current_unbonded_periodic_atoms -= remove

    def _judge_region(self, x1, y1, z1, x2, y2, z2, region, mol, basis_mol, res_len, box):
        if (x2, y2, z2) not in region:
            return
        mol |= basis_mol
        for residue in mol.residues[res_len:]:
            for atom in residue.atoms:
                atom.x = atom.x - x1 + x2
                atom.y = atom.y - y1 + y2
                atom.z = atom.z - z1 + z2
        self._process_periodic_bonds(mol, mol.residues[-1], box)

    def create(self, box, region, mol=None):
        if not isinstance(box, (BlockRegion, PrismRegion)) or box.side == "out":
            raise ValueError("Box should only be a BlockRegion or PrismRegion with side == 'in' !")
        if mol is None:
            mol = Molecule("unnamed")
        if isinstance(box, PrismRegion):
            box_length, box_angle = _length_angle_from_basis_vectors(box.l1, box.l2, box.l3)
            mol.box_length = box_length
            mol.box_angle = box_angle
            vertices = np.array([
                box.l0, box.l0 + box.l1, box.l0 + box.l2, box.l0 + box.l3,
                box.l0 + box.l1 + box.l2, box.l0 + box.l1 + box.l3,
                box.l0 + box.l2 + box.l3, box.l0 + box.l1 + box.l2 + box.l3,
            ])
            x_low, y_low, z_low = np.min(vertices, axis=0)
            x_high, y_high, z_high = np.max(vertices, axis=0)
            periodic_box = BlockRegion(x_low, y_low, z_low, x_high, y_high, z_high, boundary=True)
        else:
            mol.box_length = [box.x_high - box.x_low, box.y_high - box.y_low, box.z_high - box.z_low]
            x_low, y_low, z_low = box.x_low, box.y_low, box.z_low
            x_high, y_high, z_high = box.x_high, box.y_high, box.z_high
            periodic_box = box
        basis_mol = self.basis_molecule
        res_len = -1
        if isinstance(basis_mol, Molecule):
            res_len = -len(basis_mol.residues)
        basis_vectors = np.array(_basis_vectors_from_length_and_angle(
            self.cell_length[0] + self.spacing[0],
            self.cell_length[1] + self.spacing[1],
            self.cell_length[2] + self.spacing[2],
            self.cell_angle[0],
            self.cell_angle[1],
            self.cell_angle[2],
        )) * self.scale
        basis_positions = np.array([
            [
                basis_vectors[0][0] * basis[0] + basis_vectors[1][0] * basis[1] + basis_vectors[2][0] * basis[2],
                basis_vectors[1][1] * basis[1] + basis_vectors[2][1] * basis[2],
                basis_vectors[2][2] * basis[2],
            ]
            for basis in self.basis_position
        ])
        x_init = x_low + self.origin[0]
        y_init = y_low + self.origin[1]
        z0 = z_low + self.origin[2]
        x1, y1, z1 = np.min([[atom.x, atom.y, atom.z] for atom in basis_mol.atoms], axis=0)
        while z0 < z_high:
            y0 = y_init
            while y0 < y_high:
                x0 = x_init
                while x0 < x_high:
                    for basis in basis_positions:
                        x2 = basis[0] + x0
                        y2 = basis[1] + y0
                        z2 = basis[2] + z0
                        self._judge_region(x1, y1, z1, x2, y2, z2, region, mol, basis_mol, res_len, periodic_box)
                    x0 += basis_vectors[0][0]
                x_init += basis_vectors[1][0]
                x_init %= basis_vectors[0][0]
                y0 += basis_vectors[1][1]
            x_init += basis_vectors[2][0]
            x_init %= basis_vectors[0][0]
            y_init += basis_vectors[2][1]
            y_init %= basis_vectors[1][1]
            z0 += basis_vectors[2][2]
        return mol


SIMPLE_CUBIC_LATTICE = Lattice("template:sc", basis_position=[[0, 0, 0]])
BODY_CENTERED_CUBIC_LATTICE = Lattice("template:bcc", basis_position=[[0, 0, 0], [0.5, 0.5, 0.5]])
FACE_CENTERED_CUBIC_LATTICE = Lattice("template:fcc", basis_position=[[0, 0, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0]])
HEXAGONAL_CLOSE_PACKED_LATTICE = Lattice("template:hcp", basis_position=[[0, 0, 0], [1 / 3, 1 / 3, 0.5]], cell_angle=[90, 90, 60], cell_length=[1, 1, 2 / 3 * math.sqrt(6)])
DIAMOND_LATTICE = Lattice("template:diamond", basis_position=[[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0], [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25]])

Impose_Bond = impose_bond
Impose_Angle = impose_angle
Impose_Dihedral = impose_dihedral
Main_Axis_Rotate = main_axis_rotate
Get_Peptide_From_Sequence = get_peptide_from_sequence
H_Mass_Repartition = h_mass_repartition
Solvent_Replace = solvent_replace
Sort_Atoms_By = sort_atoms_by
Optimize = optimize
