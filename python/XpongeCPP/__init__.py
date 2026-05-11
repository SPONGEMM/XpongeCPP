"""Python compatibility layer for the XpongeCPP C++ core."""

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
    load_frcmod,
    load_charmm_parameter_file,
    load_charmm_topology_file,
    load_edip_parameter_file,
    load_gromacs_topology_file,
    load_mol2,
    load_opls_itp_file,
    load_parmdat,
    load_pdb,
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
    save_mol2,
    save_sponge_input,
    set_box_padding,
    template_atom_count,
)

__version__ = "0.1.0"


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


Load_PDB = load_pdb
LoadPDB = load_pdb
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


def get_assignment_from_pubchem(name, total_charge=None, **kwargs):
    try:
        import pubchempy as pcp
    except ImportError as exc:
        raise ImportError("PubChemPy is required for get_assignment_from_pubchem") from exc

    compounds = pcp.get_compounds(name, "name")
    if not compounds:
        raise ValueError(f"PubChem query returned no compounds: {name}")
    smiles = compounds[0].isomeric_smiles or compounds[0].canonical_smiles
    if not smiles:
        raise ValueError(f"PubChem compound has no SMILES: {name}")
    assignment = get_assignment_from_smiles(smiles, total_charge=total_charge, **kwargs)
    assignment.name = str(getattr(compounds[0], "cid", "PUBCHEM"))
    return assignment


def get_assignment_from_cif(source):
    text = source.read() if hasattr(source, "read") else open(source, encoding="utf-8").read()
    rows = []
    headers = []
    in_loop = False
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line == "loop_":
            in_loop = True
            headers = []
            continue
        if in_loop and line.startswith("_atom_site."):
            headers.append(line)
            continue
        if in_loop and headers:
            parts = line.split()
            if len(parts) >= len(headers):
                rows.append(dict(zip(headers, parts)))
            elif line.startswith("_"):
                in_loop = False
    if not rows:
        raise ValueError("CIF input does not contain a simple _atom_site loop")

    assignment = Assign("CIF")
    for row in rows:
        element = row.get("_atom_site.type_symbol", "").strip("'\"")
        name = row.get("_atom_site.label_atom_id", element).strip("'\"")
        x = float(row.get("_atom_site.Cartn_x", "0").strip("'\""))
        y = float(row.get("_atom_site.Cartn_y", "0").strip("'\""))
        z = float(row.get("_atom_site.Cartn_z", "0").strip("'\""))
        assignment.add_atom(element, x, y, z, name)
    assignment.determine_connectivity(1.2)
    return assignment


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
    "Save_Mol2",
    "Load_Gromacs_Topology_File",
    "LoadGromacsTopologyFile",
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
