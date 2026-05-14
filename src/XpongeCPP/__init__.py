"""Python compatibility layer for the XpongeCPP C++ core."""

import math
import re
from io import StringIO
from importlib.resources import files
from collections import OrderedDict
from pathlib import Path

import numpy as np

from ._core import (
    Atom,
    Assign,
    Molecule,
    Residue,
    ResidueType,
    add_ions,
    add_molecule,
    add_solvent_box as _core_add_solvent_box,
    get_assignment_from_mol2,
    get_assignment_from_pdb,
    get_assignment_from_residuetype,
    get_assignment_from_xyz,
    get_template_molecule,
    has_template,
    implemented_gaff_assign_types,
    load_coordinate,
    load_frcmod,
    load_gro,
    load_charmm_parameter_file,
    load_charmm_topology_file,
    load_edip_parameter_file,
    load_gromacs_topology_file,
    load_molpsf,
    load_mol2 as _core_load_mol2,
    load_opls_itp_file,
    load_parmdat,
    load_pdb,
    load_rst7,
    load_sw_parameter_file,
    merge_dual_topology,
    merge_force_field,
    molecule_from_residuetype,
    reorder_atoms_by_template,
    register_ff14sb,
    register_amber_frcmod_file,
    register_amber_bond_parameter,
    register_amber_lj_parameter,
    register_amber_cmap_parameter,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_text,
    register_residue_templates_from_mol2_file,
    register_residue_template_alias,
    register_pdb_residue_alias_mapping,
    register_pdb_residue_name_mapping,
    register_his_mapping,
    register_template_molecule_from_mol2_file,
    register_template_virtual_atom2,
    register_tip3p,
    registered_template_names,
    save_pdb,
    save_gro,
    save_mol2,
    save_sponge_input,
    set_box_padding,
    configure_residue_template_connect_atom,
    configure_residue_template_head,
    configure_residue_template_tail,
    replace_residues,
    template_atom_count,
)
from .assign import AssignRule
from .gromacs import GlobalSetting, GromacsTopologyIterator, load_ffitp, load_molitp

__version__ = "0.1.0"

_core_determine_connectivity = Assign.determine_connectivity
_core_determine_bond_order = Assign.determine_bond_order
_core_determine_atom_type = Assign.determine_atom_type
_core_save_as_mol2 = Assign.save_as_mol2
_core_save_as_pdb = Assign.save_as_pdb

_CONNECTIVITY_RADII = {
    "H": 0.35,
    "C": 0.73,
    "N": 0.66,
    "O": 0.69,
    "F": 0.68,
    "P": 1.06,
    "S": 1.02,
    "Cl": 0.99,
    "Br": 1.14,
    "I": 1.33,
}

_legacy_template_metadata = {}


def _tpacm4_tables():
    base = files("XpongeCPP.data.assign.tpacm4")
    return (base.joinpath("ATOMTYPE.dat").read_text(), base.joinpath("CHARGE.dat").read_text())


def _assignment_to_rdkit_mol(assignment):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("RDKit is required for Gasteiger charge calculation") from exc

    mol = Chem.RWMol()
    for atom_index, element in enumerate(assignment.atoms):
        atom = Chem.Atom(element)
        formal_charges = getattr(assignment, "formal_charges", [])
        if atom_index < len(formal_charges):
            atom.SetFormalCharge(int(formal_charges[atom_index]))
        mol.AddAtom(atom)
    bond_type = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for atom1, bonded_atoms in enumerate(assignment.bonds):
        for atom2, order in bonded_atoms.items():
            if atom1 < atom2:
                mol.AddBond(atom1, int(atom2), bond_type.get(int(order), Chem.BondType.UNSPECIFIED))
    rdmol = mol.GetMol()
    conf = Chem.Conformer(assignment.atom_count)
    for atom in range(assignment.atom_count):
        x, y, z = assignment.coordinates[atom]
        conf.SetAtomPosition(atom, (x, y, z))
    rdmol.AddConformer(conf)
    flags = Chem.rdmolops.SanitizeFlags
    Chem.SanitizeMol(rdmol, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
    return rdmol


def _rdmol_to_assignment(rdmol):
    assignment = Assign()
    conf = rdmol.GetConformer()
    for atom in rdmol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        assignment.add_atom(atom.GetSymbol(), pos.x, pos.y, pos.z)
        if atom.GetFormalCharge():
            assignment.set_formal_charge(assignment.atom_count - 1, atom.GetFormalCharge())
    has_unknown_bond = False
    from rdkit import Chem
    for bond in rdmol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.UNSPECIFIED:
            order = -1
        elif bond_type == Chem.BondType.SINGLE:
            order = 1
        elif bond_type == Chem.BondType.DOUBLE:
            order = 2
        elif bond_type == Chem.BondType.TRIPLE:
            order = 3
        elif bond_type == Chem.BondType.AROMATIC:
            order = -1
            has_unknown_bond = True
        else:
            raise NotImplementedError(f"Unknown RDKit bond type {bond_type}")
        assignment.add_bond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order)
    assignment.determine_ring_and_bond_type()
    if has_unknown_bond:
        assignment.determine_bond_order()
    return assignment


def _assign_calculate_charge(self, method, **parameters):
    method = method.upper()
    if method == "TPACM4":
        atom_type_table, charge_table = _tpacm4_tables()
        charge = parameters.get("charge", int(round(sum(self.formal_charges))))
        self._calculate_tpacm4(atom_type_table, charge_table, int(charge))
        return None
    if method == "GASTEIGER":
        try:
            from rdkit.Chem import rdPartialCharges
        except ImportError as exc:
            raise ImportError("RDKit is required for Gasteiger charge calculation") from exc
        rdmol = _assignment_to_rdkit_mol(self)
        rdPartialCharges.ComputeGasteigerCharges(rdmol)
        self.set_charges([float(atom.GetProp("_GasteigerCharge")) for atom in rdmol.GetAtoms()])
        return None
    if method == "RESP":
        from .assign import resp
        self.set_charges(
            resp.resp_fit(
                self,
                basis=parameters.get("basis", "6-31g*"),
                opt=parameters.get("opt", False),
                charge=parameters.get("charge", int(round(sum(self.formal_charges)))),
                spin=parameters.get("spin", 0),
                extra_equivalence=parameters.get("extra_equivalence", []),
                grid_density=parameters.get("grid_density", 6),
                grid_cell_layer=parameters.get("grid_cell_layer", 4),
                a1=parameters.get("a1", 0.0005),
                a2=parameters.get("a2", 0.001),
                two_stage=parameters.get("two_stage", True),
                only_esp=parameters.get("only_esp", False),
                radius=parameters.get("radius", None),
            )
        )
        return None
    raise ValueError("methods should be one of the following: 'RESP', 'GASTEIGER', 'TPACM4' (case-insensitive)")


def _first_neighbor(assignment, atom):
    return next(iter(assignment.bonds[atom]))


def _is_carboxylic_acid_hydrogen(assignment, atom):
    if assignment.atoms[atom] != "H":
        return False
    oxygen = _first_neighbor(assignment, atom)
    if assignment.atoms[oxygen] != "O":
        return False
    carbon = None
    for neighbor in assignment.bonds[oxygen]:
        if neighbor != atom:
            carbon = neighbor
            break
    if carbon is None or not assignment.atoms[carbon] == "C" or len(assignment.bonds[carbon]) != 3:
        return False
    for neighbor, order in assignment.bonds[carbon].items():
        if neighbor != oxygen and assignment.atoms[neighbor] == "O" and int(order) == 2:
            return True
    return False


def _is_carboxylate_oxygen(assignment, atom):
    if assignment.atoms[atom] != "O" or assignment.formal_charges[atom] != -1:
        return False
    carbon = _first_neighbor(assignment, atom)
    for neighbor, order in assignment.bonds[carbon].items():
        if neighbor != atom and assignment.atoms[neighbor] == "O" and int(order) == 2:
            return True
    return False


def _is_phenol_hydrogen(assignment, atom):
    if assignment.atoms[atom] != "H":
        return False
    oxygen = _first_neighbor(assignment, atom)
    if assignment.atoms[oxygen] != "O":
        return False
    for neighbor in assignment.bonds[oxygen]:
        if neighbor != atom and assignment.has_atom_marker(neighbor, "AR0"):
            return True
    return False


def _is_phenolate_oxygen(assignment, atom):
    if assignment.atoms[atom] != "O" or assignment.formal_charges[atom] != -1:
        return False
    return assignment.has_atom_marker(_first_neighbor(assignment, atom), "AR0")


def _hydrogen_position_for(assignment, atom):
    x = y = z = 0.0
    coord = assignment.coordinates[atom]
    count = 0
    for neighbor in assignment.bonds[atom]:
        ncoord = assignment.coordinates[neighbor]
        dx = coord[0] - ncoord[0]
        dy = coord[1] - ncoord[1]
        dz = coord[2] - ncoord[2]
        norm = (dx * dx + dy * dy + dz * dz) ** 0.5 or 1.0
        x += coord[0] + dx / norm
        y += coord[1] + dy / norm
        z += coord[2] + dz / norm
        count += 1
    if count == 0:
        return coord
    return x / count, y / count, z / count


def _assign_set_ph(self, ph):
    self.kekulize()
    to_delete = []
    to_add = []
    for atom in range(self.atom_count):
        if _is_carboxylic_acid_hydrogen(self, atom) and 4.0 < ph:
            to_delete.append(atom)
        elif _is_phenol_hydrogen(self, atom) and 10.0 < ph:
            to_delete.append(atom)
        elif self.atoms[atom] == "H" and self.atoms[_first_neighbor(self, atom)] == "O" and 15.9 < ph:
            to_delete.append(atom)
        elif _is_carboxylate_oxygen(self, atom) and 4.0 > ph:
            to_add.append(atom)
        elif _is_phenolate_oxygen(self, atom) and 10.0 > ph:
            to_add.append(atom)
        elif self.atoms[atom] == "O" and self.formal_charges[atom] == -1 and 15.9 > ph:
            to_add.append(atom)

    for atom in to_add:
        x, y, z = _hydrogen_position_for(self, atom)
        self.add_atom("H", x, y, z)
        self.add_bond(self.atom_count - 1, atom, 1)
    for atom in sorted(set(to_delete), reverse=True):
        self.delete_atom(atom)
    self.determine_bond_order(True, None)
    return int(round(sum(self.formal_charges)))


def _distance(coord1, coord2):
    dx = coord1[0] - coord2[0]
    dy = coord1[1] - coord2[1]
    dz = coord1[2] - coord2[2]
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _assign_determine_connectivity(self, simple_cutoff=None, tolerance=1.0):
    if simple_cutoff is not None:
        cutoff = float(simple_cutoff)
        for i, coord_i in enumerate(self.coordinates):
            for j in range(i + 1, self.atom_count):
                if _distance(coord_i, self.coordinates[j]) < cutoff:
                    self.add_bond(i, j, -1)
        return None
    for i, coord_i in enumerate(self.coordinates):
        for j in range(i + 1, self.atom_count):
            distance = _distance(coord_i, self.coordinates[j])
            if distance == 0:
                continue
            rij = _CONNECTIVITY_RADII.get(self.atoms[i], 1.25) + _CONNECTIVITY_RADII.get(self.atoms[j], 1.25)
            if distance <= 1.5:
                factor = 1 - 0.15
            elif distance <= 1.9:
                factor = 1 - 0.11
            elif distance <= 2.05:
                factor = 1 - 0.09
            else:
                factor = 1 - 0.08
            factor /= tolerance
            ratio = rij / distance
            if factor < ratio < 2:
                self.add_bond(i, j, -1)
    return None


def _assign_determine_bond_order(self, max_step=2000, max_stat=20000, penalty_scores=None,
                                 check_formal_charge=True, total_charge=None, extra_criteria=None):
    if isinstance(max_step, bool):
        check_formal_charge = max_step
        total_charge = max_stat
        max_step = 2000
        max_stat = 20000
    if penalty_scores is None and extra_criteria is None and max_step == 2000 and max_stat == 20000:
        return _core_determine_bond_order(self, check_formal_charge, total_charge)
    penalties = _normalise_penalty_scores(self, penalty_scores)
    return self._determine_bond_order_custom(
        check_formal_charge,
        total_charge,
        max_step,
        max_stat,
        penalties,
        extra_criteria,
    )


def _normalise_penalty_scores(assign, penalty_scores):
    if penalty_scores is None:
        return []
    if len(penalty_scores) != assign.atom_count:
        raise ValueError("penalty_scores should have one entry per atom")
    result = []
    for entry in penalty_scores:
        if isinstance(entry, dict):
            items = list(entry.items())
        else:
            items = list(entry)
        if not items:
            raise ValueError("penalty_scores entries should not be empty")
        result.append([(int(valence), int(penalty)) for valence, penalty in items])
    return result


def _phmodel_type(self, atom):
    if _is_phenolate_oxygen(self, atom):
        return "B-phenol"
    if _is_phenol_hydrogen(self, atom):
        return "A-phenol"
    if _is_carboxylate_oxygen(self, atom):
        return "B-carboxylic"
    if _is_carboxylic_acid_hydrogen(self, atom):
        return "A-carboxylic"
    if self.atoms[atom] == "O" and self.formal_charges[atom] == -1:
        return "B-alcohol"
    if self.atoms[atom] == "H" and self.bonds[atom] and self.atoms[_first_neighbor(self, atom)] == "O":
        return "A-alcohol"
    return "N"


def _assign_determine_atom_type(self, rule):
    if str(rule).lower() == "phmodel":
        self.kekulize()
        return [_phmodel_type(self, atom) for atom in range(self.atom_count)]
    if isinstance(rule, AssignRule):
        return _assign_determine_atom_type_from_rule(self, rule)
    if str(rule).lower() not in {"gaff", "gaff2", "sybyl"} and str(rule) in AssignRule.all:
        return _assign_determine_atom_type_from_rule(self, AssignRule.all[str(rule)])
    return _core_determine_atom_type(self, rule)


def _assign_determine_atom_type_from_rule(assign, rule):
    if not rule.built:
        rule.rules = OrderedDict(sorted(rule.rules.items(), key=lambda item: rule.priority[item[0]]))
        rule.built = True
    backup = list(assign.atom_types)
    if rule.pre_action:
        rule.pre_action(assign)
    assigned_types = []
    for atom in range(assign.atom_count):
        for atom_type, judge in rule.rules.items():
            if judge(atom, assign):
                assigned_types.append(atom_type)
                assign.set_atom_type(atom, atom_type)
                break
        else:
            raise KeyError(f"No atom type found for assignment {assign.name} of atom #{atom}")
    if rule.post_action:
        rule.post_action(assign)
    if rule.pure_string:
        for atom, atom_type in enumerate(backup):
            assign.set_atom_type(atom, atom_type)
        return assigned_types
    return None


def _assign_determine_equal_atoms(self):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("RDKit is required for determine_equal_atoms") from exc
    smiles = []
    for atom in range(self.atom_count):
        mol = _assignment_to_rdkit_mol(self)
        mol.GetAtoms()[atom].SetIsotope(1)
        smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
    group = {i: i for i in range(self.atom_count)}
    for i in range(self.atom_count):
        if group[i] == i:
            for j in range(i + 1, self.atom_count):
                if smiles[j] == smiles[i]:
                    group[j] = i
    out = []
    realmap = {}
    for atom in range(self.atom_count):
        if group[atom] == atom:
            realmap[atom] = len(out)
            out.append([atom])
        else:
            out[realmap[group[atom]]].append(atom)
    return [atoms for atoms in out if len(atoms) > 1]


def _assign_uff_optimize(self):
    try:
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for uff_optimize") from exc
    mol = _assignment_to_rdkit_mol(self)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for atom in range(self.atom_count):
        pos = conf.GetAtomPosition(atom)
        self.set_coordinate(atom, pos.x, pos.y, pos.z)


def _assign_uff_energy(self):
    try:
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for uff_energy") from exc
    mol = _assignment_to_rdkit_mol(self)
    ff = AllChem.UFFGetMoleculeForceField(mol)
    return ff.CalcEnergy()


def _assign_save_as_mol2(self, filename, residue_name="MOL", atomtype="sybyl"):
    if atomtype and str(atomtype).lower() != "sybyl":
        raise ValueError(f"save_as_mol2 only supports atomtype='sybyl'; got {atomtype!r}")
    if atomtype:
        self.determine_atom_type("sybyl")
    return _core_save_as_mol2(self, filename, residue_name)


Assign.calculate_charge = _assign_calculate_charge
Assign.Calculate_Charge = _assign_calculate_charge
Assign.set_ph = _assign_set_ph
Assign.Set_PH = _assign_set_ph
Assign.determine_connectivity = _assign_determine_connectivity
Assign.Determine_Connectivity = _assign_determine_connectivity
Assign.determine_bond_order = _assign_determine_bond_order
Assign.Determine_Bond_Order = _assign_determine_bond_order
Assign.determine_atom_type = _assign_determine_atom_type
Assign.Determine_Atom_Type = _assign_determine_atom_type
Assign.determine_equal_atoms = _assign_determine_equal_atoms
Assign.Determine_Equal_Atoms = _assign_determine_equal_atoms
Assign.uff_optimize = _assign_uff_optimize
Assign.UFF_Optimize = _assign_uff_optimize
Assign.uff_energy = property(_assign_uff_energy)
Assign.save_as_mol2 = _assign_save_as_mol2
Assign.Save_As_Mol2 = _assign_save_as_mol2


def Add_Solvent_Box(molecule, solvent, distance, tolerance=2.5, n_solvent=None, seed=0):
    solvent_molecule = solvent if isinstance(solvent, Molecule) else _single_residue_molecule(solvent, "solvent")
    return _core_add_solvent_box(
        molecule,
        solvent_molecule,
        distance,
        tolerance,
        n_solvent,
        seed,
    )


def Add_Ions(molecule, counts, seed=0, solvent="WAT"):
    return add_ions(molecule, counts, seed, solvent)


def Add_Molecule(molecule, other):
    return add_molecule(molecule, other)


def Set_Box_Padding(molecule, padding=0.5, center=True):
    return set_box_padding(molecule, padding, center)


def Save_SPONGE_Input(molecule, prefix=None, dirname="."):
    return save_sponge_input(molecule, "" if prefix is None else prefix, dirname)


def Save_PDB(molecule, filename, write_cryst1=True):
    return save_pdb(molecule, filename, write_cryst1)


def Save_Mol2(molecule, filename):
    return save_mol2(molecule, filename)


def Save_GRO(molecule, filename):
    return save_gro(molecule, filename)


add_solvent_box = Add_Solvent_Box
Save_Sponge_Input = Save_SPONGE_Input


class _LegacyResidueTypeHandle:
    def __init__(self, name):
        if not has_template(name):
            raise KeyError(f"ResidueType {name!r} is not registered")
        self._name = str(name)

    @property
    def name(self):
        return self._name

    @property
    def head(self):
        return _legacy_template_metadata.get(self._name, {}).get("head")

    @head.setter
    def head(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["head"] = value
        if value:
            configure_residue_template_head(self._name, str(value))

    @property
    def tail(self):
        return _legacy_template_metadata.get(self._name, {}).get("tail")

    @tail.setter
    def tail(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["tail"] = value
        if value:
            configure_residue_template_tail(self._name, str(value))

    def __repr__(self):
        return f"<LegacyResidueTypeHandle {self._name}>"


def _legacy_get_residuetype(name):
    if not has_template(name):
        raise KeyError(f"ResidueType {name!r} is not registered")
    return _LegacyResidueTypeHandle(name)


def load_mol2(source, as_template=False):
    if not as_template:
        return _core_load_mol2(source)
    if isinstance(source, (str, Path)):
        register_residue_templates_from_mol2_file(str(source))
        return _core_load_mol2(str(source))
    if hasattr(source, "read"):
        text = source.read()
        register_residue_templates_from_mol2_text(text)
        return _core_load_mol2(StringIO(text))
    raise TypeError("load_mol2(..., as_template=True) expects a path-like object or a readable text stream")


def _coerce_atom_index(atom):
    if isinstance(atom, (int, np.integer)):
        return int(atom)
    if hasattr(atom, "index"):
        return int(atom.index)
    raise TypeError("atom references should be atom objects or integer indices")


class _AtomIndexProxy:
    def __getitem__(self, atom):
        return _coerce_atom_index(atom)

    def get(self, atom, default=None):
        try:
            return _coerce_atom_index(atom)
        except Exception:
            return default


_core_molecule_add_residue_link = Molecule.add_residue_link


def _legacy_add_residue_link(self, atom1, atom2):
    return _core_molecule_add_residue_link(self, _coerce_atom_index(atom1), _coerce_atom_index(atom2))


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


ResidueType.get_type = staticmethod(_legacy_get_residuetype)
ResidueType.Get_Type = staticmethod(_legacy_get_residuetype)
Molecule.add_residue_link = _legacy_add_residue_link
Molecule.Add_Residue_Link = _legacy_add_residue_link
Molecule.atom_index = property(lambda self: _AtomIndexProxy())
Molecule.add_missing_atoms = _legacy_add_missing_atoms
Molecule.Add_Missing_Atoms = _legacy_add_missing_atoms
Residue.unterminal = lambda self: self
Residue.Unterminal = Residue.unterminal
Residue.UnTerminal = Residue.unterminal


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
    import os
    import shutil
    import subprocess
    import tempfile

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
            Save_SPONGE_Input(mol, prefix=prefix, dirname=tempdir)
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
    v1 = np.asarray(l1, dtype=float)
    v2 = np.asarray(l2, dtype=float)
    v3 = np.asarray(l3, dtype=float)
    a = np.linalg.norm(v1)
    b = np.linalg.norm(v2)
    c = np.linalg.norm(v3)
    alpha = math.degrees(math.acos(np.dot(v2, v3) / b / c))
    beta = math.degrees(math.acos(np.dot(v1, v3) / a / c))
    gamma = math.degrees(math.acos(np.dot(v1, v2) / a / b))
    return [a, b, c], [alpha, beta, gamma]


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


Load_PDB = load_pdb
LoadPDB = load_pdb
Load_Coordinate = load_coordinate
LoadCoordinate = load_coordinate
Load_RST7 = load_rst7
LoadRST7 = load_rst7
Load_GRO = load_gro
LoadGRO = load_gro
Load_MolPSF = load_molpsf
Load_FFITP = load_ffitp
LoadFFITP = load_ffitp
Load_MolITP = load_molitp
LoadMolITP = load_molitp
LoadMolPSF = load_molpsf
Load_Mol2 = load_mol2
LoadMOL2 = load_mol2
Load_Gromacs_Topology_File = load_gromacs_topology_file
LoadGromacsTopologyFile = load_gromacs_topology_file
Load_OPLS_ITP_File = load_opls_itp_file
LoadOPLSITPFile = load_opls_itp_file
Load_CHARMM_Parameter_File = load_charmm_parameter_file
LoadCHARMMParameterFile = load_charmm_parameter_file
Load_CHARMM_Topology_File = load_charmm_topology_file
LoadCHARMMTopologyFile = load_charmm_topology_file
Load_SW_Parameter_File = load_sw_parameter_file
LoadSWParameterFile = load_sw_parameter_file
Load_EDIP_Parameter_File = load_edip_parameter_file
LoadEDIPParameterFile = load_edip_parameter_file
Get_Assignment_From_Mol2 = get_assignment_from_mol2
Get_Assignment_From_XYZ = get_assignment_from_xyz
Get_Assignment_From_PDB = get_assignment_from_pdb
Get_Assignment_From_ResidueType = get_assignment_from_residuetype
GetAssignmentFromMol2 = get_assignment_from_mol2
GetAssignmentFromXYZ = get_assignment_from_xyz
GetAssignmentFromPDB = get_assignment_from_pdb
GetAssignmentFromResidueType = get_assignment_from_residuetype
Load_Frcmod = load_frcmod
Load_Parmdat = load_parmdat
AddIons = Add_Ions
AddMolecule = Add_Molecule
AddSolventBox = Add_Solvent_Box
SolventReplace = solvent_replace
SortAtomsBy = sort_atoms_by
SetBoxPadding = Set_Box_Padding
SaveSpongeInput = Save_SPONGE_Input
Save_SPONGEInput = Save_SPONGE_Input
SavePDB = Save_PDB
Save_PDB_File = Save_PDB
SaveGRO = Save_GRO
Save_GRO_File = Save_GRO
SaveMol2 = Save_Mol2


def _hy36_decode(width, field):
    text = field.strip()
    if not text:
        return None
    if text[0] in "+-" or text[0].isdigit():
        return int(text)
    digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    digits_lower = "0123456789abcdefghijklmnopqrstuvwxyz"
    lower = field[0].islower()
    digits = digits_lower if lower else digits_upper
    value = 0
    for char in field:
        value = value * 36 + digits.index(char)
    power = 36 ** (width - 1)
    if lower:
        return value + 16 * power + 10 ** width
    return value - 10 * power + 10 ** width


def _guess_element(atom_name, explicit=""):
    explicit = explicit.strip()
    if explicit:
        return explicit[0].upper() + explicit[1:].lower()
    text = atom_name.strip()
    text = "".join(char for char in text if not char.isdigit())
    if not text:
        return ""
    if len(text) >= 2 and text[:2].capitalize() in {"Cl", "Br", "Na", "Mg", "Ca", "Zn", "Fe"}:
        return text[:2].capitalize()
    return text[0].upper()


def get_assignment_from_pdb(file, only_residue="", bond_tolerance=1.0, total_charge=None):
    text = _read_text_source(file)
    assignment = Assign()
    serial_to_atom = {}
    conect = []
    has_conect = False
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if only_residue and resname != only_residue:
                continue
            serial = _hy36_decode(5, line[6:11])
            if serial is None:
                continue
            atom_name = line[12:16].strip()
            element = _guess_element(atom_name, line[76:78] if len(line) >= 78 else "")
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            assignment.name = resname
            serial_to_atom[serial] = assignment.atom_count
            assignment.add_atom(element, x, y, z, atom_name)
        elif line.startswith("CONECT"):
            has_conect = True
            atom = _hy36_decode(5, line[6:11])
            if atom is None:
                continue
            for start in range(11, min(len(line), 31), 5):
                field = line[start:start + 5]
                if not field.strip():
                    continue
                bonded_atom = _hy36_decode(5, field)
                if bonded_atom is not None:
                    conect.append((atom, bonded_atom))
    if assignment.atom_count == 0:
        raise OSError("The input is not a pdb file")
    for atom, bonded_atom in conect:
        if atom in serial_to_atom and bonded_atom in serial_to_atom:
            atom1 = serial_to_atom[atom]
            atom2 = serial_to_atom[bonded_atom]
            if atom1 < atom2:
                assignment.add_bond(atom1, atom2, 1)
    if not has_conect:
        assignment.determine_connectivity(tolerance=bond_tolerance)
    assignment.determine_bond_order(total_charge=total_charge)
    return assignment


def get_assignment_from_xyz(file, bond_tolerance=1.0, total_charge=None):
    text = _read_text_source(file)
    lines = text.splitlines()
    if not lines:
        raise OSError("The input is not a xyz file")
    atom_numbers = int(lines[0].strip())
    assignment = Assign()
    assignment.name = lines[1].strip()
    for index, line in enumerate(lines[2:2 + atom_numbers]):
        atom_name, x, y, z = line.split()[:4]
        assignment.add_atom(atom_name, float(x), float(y), float(z), f"{atom_name}{index + 1}")
    assignment.determine_connectivity(tolerance=bond_tolerance)
    assignment.determine_bond_order(total_charge=total_charge)
    return assignment


def get_assignment_from_smiles(smiles):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for get_assignment_from_smiles") from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return _rdmol_to_assignment(mol)


def get_assignment_from_pubchem(parameter, keyword="name", total_charge=None, **kwargs):
    try:
        import pubchempy as pcp
    except ImportError as exc:
        raise ImportError("PubChemPy is required for get_assignment_from_pubchem") from exc

    compounds = pcp.get_compounds(parameter, keyword, record_type="3d")
    if not compounds:
        try:
            raise pcp.NotFoundError
        except AttributeError as exc:
            raise ValueError(f"PubChem query returned no compounds: {parameter}") from exc
    if len(compounds) != 1:
        raise NotImplementedError
    compound = compounds[0]
    assignment = Assign(str(getattr(compound, "cid", "PUBCHEM")))
    for atom in compound.atoms:
        assignment.add_atom(atom.element, atom.x, atom.y, atom.z)
    for bond in compound.bonds:
        assignment.add_bond(bond.aid1 - 1, bond.aid2 - 1, bond.order)
    assignment.determine_ring_and_bond_type()
    return assignment


def _read_text_source(source):
    if hasattr(source, "read"):
        return source.read()
    if isinstance(source, str) and ("\n" in source or source.lstrip().startswith("data_")):
        return source
    with open(source, encoding="utf-8") as handle:
        return handle.read()


def _cif_float(value):
    value = value.strip("'\"")
    if "(" in value:
        value = value.split("(", 1)[0]
    return float(value)


def _cif_loops(text):
    loops = []
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        headers = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1
        rows = []
        while index < len(lines):
            line = lines[index].strip()
            if not line or line.startswith("#"):
                index += 1
                continue
            if line == "loop_" or line.startswith("_") or line.startswith("data_"):
                break
            parts = line.split()
            if len(parts) >= len(headers):
                rows.append(dict(zip(headers, parts)))
            index += 1
        loops.append((headers, rows))
    return loops


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


def _parse_cif_symops(text, lattice_info):
    match = re.search(r"(_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz)\s+(.+?)(?!_)\n(_\S+|loop_\S*)", text, flags=re.DOTALL)
    if not match:
        return
    symops = match.group(2).replace("'", "")
    if set(symops) - set("+-,xyz0123456789\n /"):
        raise ValueError("the symmetry operator strings can only be simple math expression of x, y, z")
    symops = symops.replace("x", "1").replace("y", "1").replace("z", "1").strip()
    lattice_info["basis_position"] = [
        [eval(op) for op in symop.split(",") if op]  # pylint: disable=eval-used
        for symop in symops.split("\n")
    ]


def get_assignment_from_cif(source, total_charge=0, orthogonal_threshold=None, keep_cell_angle=True):
    text = _read_text_source(source)
    data_name = "CIF"
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("data_"):
            data_name = stripped[5:]
            break
    values = {}
    for raw in text.splitlines():
        parts = raw.strip().split(None, 1)
        if len(parts) == 2 and parts[0].startswith("_cell_"):
            values[parts[0]] = parts[1]

    lattice_info = {"scale": 1, "style": "custom"}
    basis = None
    if "_cell_length_a" in values:
        lengths = [
            _cif_float(values["_cell_length_a"]),
            _cif_float(values["_cell_length_b"]),
            _cif_float(values["_cell_length_c"]),
        ]
        lattice_info["cell_length"] = lengths
        angles = [
            _cif_float(values["_cell_angle_alpha"]),
            _cif_float(values["_cell_angle_beta"]),
            _cif_float(values["_cell_angle_gamma"]),
        ]
        if not keep_cell_angle:
            angles = [90, 90, 90]
        elif orthogonal_threshold is not None:
            angles = [90 if abs(angle - 90) < orthogonal_threshold else angle for angle in angles]
        lattice_info["cell_angle"] = angles
        basis = _basis_vectors_from_length_and_angle(*lengths, *angles)
    _parse_cif_symops(text, lattice_info)

    rows = []
    bond_rows = []
    for headers, loop_rows in _cif_loops(text):
        if "_atom_site_type_symbol" in headers:
            rows = loop_rows
        if "_geom_bond_atom_site_label_1" in headers:
            bond_rows = loop_rows
    if not rows:
        raise ValueError("CIF input does not contain a simple _atom_site loop")

    assignment = Assign(data_name)
    name_to_atom = {}
    for row in rows:
        element = row.get("_atom_site_type_symbol", row.get("_atom_site.type_symbol", "")).strip("'\"")
        name = row.get("_atom_site_label", row.get("_atom_site_label_atom_id", element)).strip("'\"")
        if "_atom_site_Cartn_x" in row or "_atom_site.Cartn_x" in row:
            x = _cif_float(row.get("_atom_site_Cartn_x", row.get("_atom_site.Cartn_x", "0")))
            y = _cif_float(row.get("_atom_site_Cartn_y", row.get("_atom_site.Cartn_y", "0")))
            z = _cif_float(row.get("_atom_site_Cartn_z", row.get("_atom_site.Cartn_z", "0")))
        elif "_atom_site_fract_x" in row:
            if basis is None:
                raise ValueError("fractional CIF coordinates require cell information")
            fx = _cif_float(row["_atom_site_fract_x"])
            fy = _cif_float(row["_atom_site_fract_y"])
            fz = _cif_float(row["_atom_site_fract_z"])
            x = fx * basis[0][0] + fy * basis[1][0] + fz * basis[2][0]
            y = fx * basis[0][1] + fy * basis[1][1] + fz * basis[2][1]
            z = fx * basis[0][2] + fy * basis[1][2] + fz * basis[2][2]
        else:
            raise ValueError("There is no atom position found in CIF input")
        name_to_atom[name] = assignment.atom_count
        assignment.add_atom(element, x, y, z, name)
    if bond_rows:
        for row in bond_rows:
            atom1 = row["_geom_bond_atom_site_label_1"].strip("'\"")
            atom2 = row["_geom_bond_atom_site_label_2"].strip("'\"")
            if atom1 in name_to_atom and atom2 in name_to_atom:
                assignment.add_bond(name_to_atom[atom1], name_to_atom[atom2], -1)
    else:
        assignment.determine_connectivity(1.2)
    assignment.determine_bond_order(True, total_charge)
    return assignment, lattice_info


Get_Assignment_From_Smiles = get_assignment_from_smiles
Get_Assignment_From_PubChem = get_assignment_from_pubchem
Get_Assignment_From_CIF = get_assignment_from_cif
GetAssignmentFromSmiles = get_assignment_from_smiles
GetAssignmentFromPubChem = get_assignment_from_pubchem
GetAssignmentFromCIF = get_assignment_from_cif
Get_Assignment_From_XYZ = get_assignment_from_xyz
Get_Assignment_From_PDB = get_assignment_from_pdb
GetAssignmentFromXYZ = get_assignment_from_xyz
GetAssignmentFromPDB = get_assignment_from_pdb

__all__ = [
    "Assign",
    "AssignRule",
    "GlobalSetting",
    "GromacsTopologyIterator",
    "Molecule",
    "ResidueType",
    "add_ions",
    "load_pdb",
    "load_mol2",
    "load_molpsf",
    "load_ffitp",
    "load_molitp",
    "load_coordinate",
    "load_rst7",
    "load_gro",
    "load_frcmod",
    "load_parmdat",
    "load_gromacs_topology_file",
    "load_opls_itp_file",
    "load_charmm_parameter_file",
    "load_charmm_topology_file",
    "load_sw_parameter_file",
    "load_edip_parameter_file",
    "add_solvent_box",
    "add_molecule",
    "get_assignment_from_mol2",
    "get_assignment_from_xyz",
    "get_assignment_from_pdb",
    "get_assignment_from_residuetype",
    "get_assignment_from_cif",
    "get_assignment_from_smiles",
    "get_assignment_from_pubchem",
    "set_box_padding",
    "save_sponge_input",
    "save_pdb",
    "save_gro",
    "save_mol2",
    "register_ff14sb",
    "register_tip3p",
    "register_amber_parmdat_file",
    "register_amber_frcmod_file",
    "register_amber_lj_parameter",
    "register_amber_cmap_parameter",
    "register_amber_bond_parameter",
    "register_residue_templates_from_mol2_text",
    "register_residue_templates_from_mol2_file",
    "register_residue_template_alias",
    "register_pdb_residue_alias_mapping",
    "register_pdb_residue_name_mapping",
    "register_his_mapping",
    "register_template_molecule_from_mol2_file",
    "register_template_virtual_atom2",
    "configure_residue_template_head",
    "configure_residue_template_tail",
    "configure_residue_template_connect_atom",
    "has_template",
    "template_atom_count",
    "registered_template_names",
    "get_template_molecule",
    "implemented_gaff_assign_types",
    "merge_dual_topology",
    "merge_force_field",
    "Add_Ions",
    "Add_Molecule",
    "Add_Solvent_Box",
    "Set_Box_Padding",
    "Save_SPONGE_Input",
    "Save_PDB",
    "Save_GRO",
    "Save_Mol2",
    "Load_Gromacs_Topology_File",
    "LoadGromacsTopologyFile",
    "Load_Coordinate",
    "LoadCoordinate",
    "Load_RST7",
    "LoadRST7",
    "Load_OPLS_ITP_File",
    "LoadOPLSITPFile",
    "Load_CHARMM_Parameter_File",
    "LoadCHARMMParameterFile",
    "Load_CHARMM_Topology_File",
    "LoadCHARMMTopologyFile",
    "Load_SW_Parameter_File",
    "LoadSWParameterFile",
    "Load_EDIP_Parameter_File",
    "LoadEDIPParameterFile",
]
