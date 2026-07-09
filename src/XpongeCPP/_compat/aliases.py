"""Registry-oriented name alias helpers for legacy API unification."""

from __future__ import annotations

LEGACY_TOP_LEVEL_ALIAS_SPECS = {
    "Impose_Bond": "impose_bond",
    "Impose_Angle": "impose_angle",
    "Impose_Dihedral": "impose_dihedral",
    "Main_Axis_Rotate": "main_axis_rotate",
    "Get_Peptide_From_Sequence": "get_peptide_from_sequence",
    "H_Mass_Repartition": "h_mass_repartition",
    "Solvent_Replace": "solvent_replace",
    "Sort_Atoms_By": "sort_atoms_by",
    "Optimize": "optimize",
    "Load_PDB": "load_pdb",
    "Load_Pdb": "load_pdb",
    "LoadPDB": "load_pdb",
    "LoadPdb": "load_pdb",
    "Load_MMCIF": "load_mmcif",
    "Load_Mmcif": "load_mmcif",
    "LoadMMCIF": "load_mmcif",
    "LoadMmcif": "load_mmcif",
    "Load_Coordinate": "load_coordinate",
    "LoadCoordinate": "load_coordinate",
    "Load_RST7": "load_rst7",
    "LoadRST7": "load_rst7",
    "Load_GRO": "load_gro",
    "Load_Gro": "load_gro",
    "LoadGRO": "load_gro",
    "LoadGro": "load_gro",
    "Load_MolPSF": "load_molpsf",
    "Load_FFITP": "load_ffitp",
    "LoadFFITP": "load_ffitp",
    "Load_MolITP": "load_molitp",
    "LoadMolITP": "load_molitp",
    "LoadMolPSF": "load_molpsf",
    "Load_Mol2": "load_mol2",
    "Load_Mol_2": "load_mol2",
    "LoadMOL2": "load_mol2",
    "LoadMol2": "load_mol2",
    "Load_Gromacs_Topology_File": "load_gromacs_topology_file",
    "LoadGromacsTopologyFile": "load_gromacs_topology_file",
    "Load_OPLS_ITP_File": "load_opls_itp_file",
    "LoadOPLSITPFile": "load_opls_itp_file",
    "Load_Parameter_From_FFITP": "load_parameter_from_ffitp",
    "LoadParameterFromFFITP": "load_parameter_from_ffitp",
    "Load_CHARMM_Parameter_File": "load_charmm_parameter_file",
    "LoadCHARMMParameterFile": "load_charmm_parameter_file",
    "Load_CHARMM_Topology_File": "load_charmm_topology_file",
    "LoadCHARMMTopologyFile": "load_charmm_topology_file",
    "Load_SW_Parameter_File": "load_sw_parameter_file",
    "LoadSWParameterFile": "load_sw_parameter_file",
    "Load_EDIP_Parameter_File": "load_edip_parameter_file",
    "LoadEDIPParameterFile": "load_edip_parameter_file",
    "Get_Assignment_From_Mol2": "get_assignment_from_mol2",
    "Get_Assignment_From_XYZ": "get_assignment_from_xyz",
    "Get_Assignment_From_PDB": "get_assignment_from_pdb",
    "Get_Assignment_From_ResidueType": "get_assignment_from_residuetype",
    "GetAssignmentFromMol2": "get_assignment_from_mol2",
    "GetAssignmentFromXYZ": "get_assignment_from_xyz",
    "GetAssignmentFromPDB": "get_assignment_from_pdb",
    "GetAssignmentFromResidueType": "get_assignment_from_residuetype",
    "Load_Frcmod": "load_frcmod",
    "Load_Parmdat": "load_parmdat",
    "AddIons": "Add_Ions",
    "AddMolecule": "Add_Molecule",
    "AddSolventBox": "Add_Solvent_Box",
    "SolventReplace": "solvent_replace",
    "SortAtomsBy": "sort_atoms_by",
    "SetBoxPadding": "Set_Box_Padding",
    "SaveSpongeInput": "Save_SPONGE_Input",
    "Save_SPONGEInput": "Save_SPONGE_Input",
    "SavePDB": "Save_PDB",
    "Save_PDB_File": "Save_PDB",
    "SaveGRO": "Save_GRO",
    "Save_GRO_File": "Save_GRO",
    "SaveMol2": "Save_Mol2",
    "Get_Assignment_From_Smiles": "get_assignment_from_smiles",
    "Get_Assignment_From_PubChem": "get_assignment_from_pubchem",
    "Get_Assignment_From_CIF": "get_assignment_from_cif",
    "GetAssignmentFromSmiles": "get_assignment_from_smiles",
    "GetAssignmentFromPubChem": "get_assignment_from_pubchem",
    "GetAssignmentFromCIF": "get_assignment_from_cif",
}


def install_aliases(namespace: dict, alias_map: dict[str, object]):
    """Install a batch of aliases into *namespace*."""

    for alias, target in alias_map.items():
        namespace[alias] = target
    return alias_map


def install_top_level_aliases(namespace: dict):
    """Install the standard legacy Xponge top-level alias set."""

    alias_map = {alias: namespace[target_name] for alias, target_name in LEGACY_TOP_LEVEL_ALIAS_SPECS.items()}
    return install_aliases(namespace, alias_map)
