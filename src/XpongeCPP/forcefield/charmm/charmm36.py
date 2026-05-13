"""Register packaged CHARMM36 parameters and residue-template semantics."""

from ... import (
    configure_residue_template_connect_atom,
    configure_residue_template_head,
    configure_residue_template_tail,
    load_gromacs_topology_file,
    register_his_mapping,
    register_pdb_residue_name_mapping,
    register_residue_template_alias,
    register_residue_templates_from_mol2_file,
)
from . import data_path


def _configure_protein_terminals():
    residues = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
    for resname in residues:
        configure_residue_template_head(resname, "N", 1.3, "CA")
        configure_residue_template_tail(resname, "C", 1.3, "CA")
        configure_residue_template_tail("N" + resname, "C", 1.3, "CA")
        configure_residue_template_head("C" + resname, "N", 1.3, "CA")
        register_pdb_residue_name_mapping("head", resname, "N" + resname)
        register_pdb_residue_name_mapping("tail", resname, "C" + resname)

    configure_residue_template_tail("ACE", "C", 1.3, "CH3")
    configure_residue_template_head("NME", "N", 1.3, "CH3")

    for source, alias in (("HIE", "HIS"), ("NHIE", "NHIS"), ("CHIE", "CHIS")):
        register_residue_template_alias(alias, source)
    register_pdb_residue_name_mapping("head", "HIS", "HIS")
    register_pdb_residue_name_mapping("tail", "HIS", "HIS")
    register_his_mapping("HIS", "HID", "HIE", "HIP")
    register_his_mapping("NHIS", "NHID", "NHIE", "NHIP")
    register_his_mapping("CHIS", "CHID", "CHIE", "CHIP")

    for resname in ("CYX", "NCYX", "CCYX"):
        configure_residue_template_connect_atom(resname, "ssbond", "SG")


def _configure_dna_terminals():
    for base in "ATCG":
        resname = "D" + base
        configure_residue_template_head(resname, "P", 1.5, "OP2")
        configure_residue_template_head(resname + "3", "P", 1.5, "OP2")
        configure_residue_template_tail(resname, "O3'", 1.5, "C3'")
        configure_residue_template_tail(resname + "5", "O3'", 1.5, "C3'")
        register_pdb_residue_name_mapping("head", resname, resname + "5")
        register_pdb_residue_name_mapping("tail", resname, resname + "3")


def _configure_rna_terminals():
    for resname in "AUCG":
        configure_residue_template_head(resname, "P", 1.5, "OP2")
        configure_residue_template_head(resname + "3", "P", 1.5, "OP2")
        configure_residue_template_tail(resname, "O3'", 1.5, "C3'")
        configure_residue_template_tail(resname + "5", "O3'", 1.5, "C3'")
        register_pdb_residue_name_mapping("head", resname, resname + "5")
        register_pdb_residue_name_mapping("tail", resname, resname + "3")


load_gromacs_topology_file(str(data_path("charmm36", "forcefield.itp")))
register_residue_templates_from_mol2_file(str(data_path("charmm36", "protein.mol2")))
register_residue_templates_from_mol2_file(str(data_path("charmm36", "dna.mol2")))
register_residue_templates_from_mol2_file(str(data_path("charmm36", "rna.mol2")))
_configure_protein_terminals()
_configure_dna_terminals()
_configure_rna_terminals()
