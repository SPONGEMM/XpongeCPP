"""Python compatibility layer for the XpongeCPP C++ core."""

from importlib.resources import files

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

__version__ = "0.1.0"


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
    bond_type = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        -1: Chem.BondType.AROMATIC,
    }
    for atom1, bonded_atoms in enumerate(assignment.bonds):
        for atom2, order in bonded_atoms.items():
            if atom1 < atom2:
                mol.AddBond(atom1, int(atom2), bond_type.get(int(order), Chem.BondType.SINGLE))
                if int(order) == -1:
                    mol.GetAtomWithIdx(atom1).SetIsAromatic(True)
                    mol.GetAtomWithIdx(int(atom2)).SetIsAromatic(True)
                    mol.GetBondBetweenAtoms(atom1, int(atom2)).SetIsAromatic(True)
    return mol.GetMol()


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
        try:
            import pyscf  # noqa: F401
        except ImportError as exc:
            raise ImportError("PySCF is required for RESP charge calculation") from exc
        raise NotImplementedError("RESP charge fitting is not yet available in XpongeCPP without the Xponge RESP backend")
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


Assign.calculate_charge = _assign_calculate_charge
Assign.Calculate_Charge = _assign_calculate_charge
Assign.set_ph = _assign_set_ph
Assign.Set_PH = _assign_set_ph


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


def get_assignment_from_smiles(smiles, total_charge=None, add_hydrogens=True, seed=20260509):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError("RDKit is required for get_assignment_from_smiles") from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    if add_hydrogens:
        mol = Chem.AddHs(mol)
    if mol.GetNumConformers() == 0:
        status = AllChem.EmbedMolecule(mol, randomSeed=int(seed))
        if status != 0:
            raise ValueError("RDKit failed to embed the SMILES molecule")
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    conf = mol.GetConformer()
    assignment = Assign("SMILES")
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        assignment.add_atom(atom.GetSymbol(), pos.x, pos.y, pos.z, atom.GetSymbol() + str(atom.GetIdx() + 1), 0.0)
    order_map = {
        Chem.BondType.SINGLE: 1,
        Chem.BondType.DOUBLE: 2,
        Chem.BondType.TRIPLE: 3,
        Chem.BondType.AROMATIC: -1,
    }
    for bond in mol.GetBonds():
        assignment.add_bond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order_map.get(bond.GetBondType(), 1))
    if total_charge is not None:
        assignment.determine_bond_order(True, total_charge)
    return assignment


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
        raise NotImplementedError("get_assignment_from_pubchem expects exactly one PubChem result")
    compound = compounds[0]
    if getattr(compound, "atoms", None) and getattr(compound, "bonds", None):
        assignment = Assign(str(getattr(compound, "cid", "PUBCHEM")))
        for atom in compound.atoms:
            assignment.add_atom(atom.element, atom.x, atom.y, atom.z)
        for bond in compound.bonds:
            assignment.add_bond(bond.aid1 - 1, bond.aid2 - 1, bond.order)
        assignment.determine_bond_order(True, total_charge)
        return assignment
    smiles = compound.isomeric_smiles or compound.canonical_smiles
    if not smiles:
        raise ValueError(f"PubChem compound has no SMILES: {parameter}")
    assignment = get_assignment_from_smiles(smiles, total_charge=total_charge, **kwargs)
    assignment.name = str(getattr(compound, "cid", "PUBCHEM"))
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
    if "_cell_length_a" in values:
        lattice_info["cell_length"] = [
            _cif_float(values["_cell_length_a"]),
            _cif_float(values["_cell_length_b"]),
            _cif_float(values["_cell_length_c"]),
        ]
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
        x = _cif_float(row.get("_atom_site_Cartn_x", row.get("_atom_site.Cartn_x", "0")))
        y = _cif_float(row.get("_atom_site_Cartn_y", row.get("_atom_site.Cartn_y", "0")))
        z = _cif_float(row.get("_atom_site_Cartn_z", row.get("_atom_site.Cartn_z", "0")))
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

__all__ = [
    "Assign",
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
