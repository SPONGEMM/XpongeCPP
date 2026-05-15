"""Python compatibility layer for the XpongeCPP C++ core."""

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
    get_assignment_from_pdb as _core_get_assignment_from_pdb,
    get_assignment_from_residuetype,
    get_assignment_from_xyz as _core_get_assignment_from_xyz,
    get_template_molecule,
    has_template,
    implemented_gaff_assign_types,
    load_coordinate,
    load_frcmod as _core_load_frcmod,
    load_gro,
    load_charmm_parameter_file,
    load_charmm_topology_file,
    load_edip_parameter_file,
    load_gromacs_topology_file,
    load_molpsf,
    load_opls_itp_file,
    load_parmdat,
    load_pdb as _core_load_pdb,
    load_rst7,
    load_sw_parameter_file,
    merge_dual_topology,
    merge_force_field,
    molecule_from_residuetype,
    reorder_atoms_by_template,
    register_ff14sb,
    clear_amber_dihedral_parameters,
    clear_amber_improper_parameters,
    register_amber_frcmod_file,
    register_amber_angle_parameter,
    register_amber_bond_parameter,
    register_amber_improper_dihedral_parameter,
    register_amber_lj_parameter,
    register_amber_nb14_scale,
    register_amber_proper_dihedral_parameter,
    register_amber_cmap_parameter,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_text as _core_register_residue_templates_from_mol2_text,
    register_residue_templates_from_mol2_file as _core_register_residue_templates_from_mol2_file,
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
    set_lj_combining_rule,
    set_box_padding,
    configure_residue_template_connect_atom,
    configure_residue_template_head,
    configure_residue_template_tail,
    replace_residues,
    template_atom_count,
)
from .assign import AssignRule
from .gromacs import GlobalSetting, GromacsTopologyIterator, load_ffitp, load_molitp
from .io_compat import (
    get_assignment_from_cif,
    get_assignment_from_pdb,
    get_assignment_from_pubchem,
    get_assignment_from_smiles,
    get_assignment_from_xyz,
)

_CoreMolecule = Molecule
_CoreResidue = Residue
_CoreResidueType = ResidueType
from ._compat.imports import (
    Generate_New_Bonded_Force_Type,
    Generate_New_Pairwise_Force_Type,
    Xdict,
    Xopen,
    Xpri,
    Xprint,
    debug,
    generate_new_bonded_force_type,
    generate_new_pairwise_force_type,
    set_global_alternative_names,
    source,
    xopen,
    xprint,
)
from ._compat.assign import install_legacy_assign_patches
from ._compat.aliases import install_top_level_aliases
from ._compat.runtime import install_legacy_runtime_patches
from ._compat.symbols import sync_template_module_globals
from ._compat.workflows import build_bonded_force, ensure_mindsponge_todo_support, get_mindsponge_system_energy
from .helper import AbstractMolecule, AtomType, Entity, ResidueLink, Type
from .helper.file import file_filter, import_python_script, pdb_filter
from .helper.math import (
    get_fibonacci_grid,
    get_rotate_matrix,
    guess_element_from_mass,
    kabsch,
)
from .process import (
    BODY_CENTERED_CUBIC_LATTICE,
    DIAMOND_LATTICE,
    FACE_CENTERED_CUBIC_LATTICE,
    HEXAGONAL_CLOSE_PACKED_LATTICE,
    SIMPLE_CUBIC_LATTICE,
    Add_Ions,
    Add_Molecule,
    Add_Solvent_Box,
    BlockRegion,
    FrustumRegion,
    IntersectRegion,
    Lattice,
    Main_Axis_Rotate,
    Optimize,
    PrismRegion,
    Region,
    Save_GRO,
    Save_Mol2,
    Save_PDB,
    Save_SPONGE_Input,
    Save_Sponge_Input,
    Set_Box_Padding,
    Solvent_Replace,
    SphereRegion,
    Sort_Atoms_By,
    UnionRegion,
    add_solvent_box,
    get_peptide_from_sequence,
    h_mass_repartition,
    impose_angle,
    impose_bond,
    impose_dihedral,
    main_axis_rotate,
    optimize,
    solvent_replace,
    sort_atoms_by,
)
from .legacy_types import _LegacyResidueTypeHandle
from .template_ops import load_mol2

__version__ = "0.1.0"
pi = np.pi
kb = 0.00198716
bar = 1.439506089041446e-5
save_mol2 = Save_Mol2
save_sponge_input = Save_SPONGE_Input


def register_residue_templates_from_mol2_file(filename):
    result = globals()["_core_register_residue_templates_from_mol2_file"](filename)
    sync_template_module_globals()
    return result


def register_residue_templates_from_mol2_text(text):
    result = globals()["_core_register_residue_templates_from_mol2_text"](text)
    sync_template_module_globals()
    return result

def load_frcmod(filename):
    set_lj_combining_rule("lorentz_berthelot")
    register_amber_nb14_scale("X", "X", 0.5, 0.833333)
    return _core_load_frcmod(filename)


def load_pdb(*args, **kwargs):
    from .forcefield import package_data_path

    set_lj_combining_rule("lorentz_berthelot")
    register_amber_nb14_scale("X", "X", 0.5, 0.833333)
    register_amber_parmdat_file(str(package_data_path("amber", "parm10.dat")))
    register_amber_frcmod_file(str(package_data_path("amber", "ff14SB.frcmod")))
    return _core_load_pdb(*args, **kwargs)


install_legacy_runtime_patches(globals())
install_legacy_assign_patches()
ensure_mindsponge_todo_support()


def load_parameter_from_ffitp(filename, folder, reset=True):
    from .forcefield.opls import load_parameter_from_ffitp as _load_parameter_from_ffitp

    return _load_parameter_from_ffitp(filename, folder, reset=reset)

install_top_level_aliases(globals())

__all__ = [
    "Assign",
    "AssignRule",
    "AbstractMolecule",
    "Atom",
    "AtomType",
    "Entity",
    "GlobalSetting",
    "Generate_New_Bonded_Force_Type",
    "Generate_New_Pairwise_Force_Type",
    "GromacsTopologyIterator",
    "Molecule",
    "Residue",
    "ResidueType",
    "ResidueLink",
    "Type",
    "Xdict",
    "Xopen",
    "Xpri",
    "Xprint",
    "xopen",
    "xprint",
    "debug",
    "file_filter",
    "set_global_alternative_names",
    "generate_new_bonded_force_type",
    "generate_new_pairwise_force_type",
    "get_fibonacci_grid",
    "get_rotate_matrix",
    "guess_element_from_mass",
    "import_python_script",
    "kabsch",
    "pdb_filter",
    "pi",
    "kb",
    "bar",
    "source",
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
    "build_bonded_force",
    "get_mindsponge_system_energy",
    "register_ff14sb",
    "register_tip3p",
    "register_amber_parmdat_file",
    "register_amber_frcmod_file",
    "register_amber_angle_parameter",
    "register_amber_lj_parameter",
    "register_amber_cmap_parameter",
    "register_amber_bond_parameter",
    "register_amber_proper_dihedral_parameter",
    "register_amber_improper_dihedral_parameter",
    "register_amber_nb14_scale",
    "clear_amber_dihedral_parameters",
    "clear_amber_improper_parameters",
    "set_lj_combining_rule",
    "load_parameter_from_ffitp",
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
    "Region",
    "UnionRegion",
    "IntersectRegion",
    "BlockRegion",
    "SphereRegion",
    "FrustumRegion",
    "PrismRegion",
    "Lattice",
    "Optimize",
    "Main_Axis_Rotate",
    "Sort_Atoms_By",
    "Solvent_Replace",
    "get_peptide_from_sequence",
    "impose_bond",
    "impose_angle",
    "impose_dihedral",
    "h_mass_repartition",
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
