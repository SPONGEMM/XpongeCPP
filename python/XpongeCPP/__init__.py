"""Python compatibility layer for the XpongeCPP C++ core."""

import math
import re
from importlib.resources import files
from collections import OrderedDict

from ._core import (
    Assign,
    Molecule,
    ResidueType,
    add_ions,
    add_molecule,
    add_solvent_box,
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
    load_mol2,
    load_opls_itp_file,
    load_parmdat,
    load_pdb,
    load_rst7,
    load_sw_parameter_file,
    merge_dual_topology,
    merge_force_field,
    register_ff14sb,
    register_amber_frcmod_file,
    register_amber_bond_parameter,
    register_amber_lj_parameter,
    register_amber_cmap_parameter,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
    register_pdb_residue_alias_mapping,
    register_pdb_residue_name_mapping,
    register_template_molecule_from_mol2_file,
    register_template_virtual_atom2,
    register_tip3p,
    save_pdb,
    save_gro,
    save_mol2,
    save_sponge_input,
    set_box_padding,
    template_atom_count,
)
from .assign import AssignRule

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
    return add_solvent_box(molecule, solvent, distance, tolerance, n_solvent, seed)


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


Load_PDB = load_pdb
LoadPDB = load_pdb
Load_Coordinate = load_coordinate
LoadCoordinate = load_coordinate
Load_RST7 = load_rst7
LoadRST7 = load_rst7
Load_GRO = load_gro
LoadGRO = load_gro
Load_MolPSF = load_molpsf
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
    "Molecule",
    "ResidueType",
    "add_ions",
    "load_pdb",
    "load_mol2",
    "load_molpsf",
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
    "register_residue_templates_from_mol2_file",
    "register_pdb_residue_alias_mapping",
    "register_pdb_residue_name_mapping",
    "register_template_molecule_from_mol2_file",
    "register_template_virtual_atom2",
    "has_template",
    "template_atom_count",
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
