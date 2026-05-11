"""Python compatibility layer for the XpongeCPP C++ core."""

from ._core import (
    Assign,
    Molecule,
    ResidueType,
    add_ions,
    add_molecule,
    add_solvent_box,
    get_assignment_from_mol2,
    get_template_molecule,
    has_template,
    implemented_gaff_assign_types,
    load_frcmod,
    load_mol2,
    load_parmdat,
    load_pdb,
    register_ff14sb,
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
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


def Add_Ions(molecule, counts, seed=0):
    return add_ions(molecule, counts, seed)


def Add_Molecule(molecule, other):
    return add_molecule(molecule, other)


def Set_Box_Padding(molecule, padding=0.5, center=True):
    return set_box_padding(molecule, padding, center)


def Save_SPONGE_Input(molecule, prefix=None, dirname="."):
    return save_sponge_input(molecule, "" if prefix is None else prefix, dirname)


def Save_PDB(molecule, filename):
    return save_pdb(molecule, filename)


def Save_Mol2(molecule, filename):
    return save_mol2(molecule, filename)


Load_PDB = load_pdb
LoadPDB = load_pdb
Load_Mol2 = load_mol2
LoadMOL2 = load_mol2
Get_Assignment_From_Mol2 = get_assignment_from_mol2
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

__all__ = [
    "Assign",
    "Molecule",
    "ResidueType",
    "add_ions",
    "load_pdb",
    "load_mol2",
    "load_frcmod",
    "load_parmdat",
    "add_solvent_box",
    "add_molecule",
    "get_assignment_from_mol2",
    "set_box_padding",
    "save_sponge_input",
    "save_pdb",
    "save_mol2",
    "register_ff14sb",
    "register_tip3p",
    "register_amber_parmdat_file",
    "register_amber_frcmod_file",
    "register_residue_templates_from_mol2_file",
    "has_template",
    "template_atom_count",
    "get_template_molecule",
    "implemented_gaff_assign_types",
    "Add_Ions",
    "Add_Molecule",
    "Add_Solvent_Box",
    "Set_Box_Padding",
    "Save_SPONGE_Input",
    "Save_PDB",
    "Save_Mol2",
]
