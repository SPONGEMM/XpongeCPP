"""Register packaged OPLS-AA/M parameters and protein residue-template semantics."""

from ... import (
    configure_residue_template_connect_atom,
    configure_residue_template_head,
    configure_residue_template_tail,
    register_his_mapping,
    register_pdb_residue_name_mapping,
    register_residue_template_alias,
    register_residue_templates_from_mol2_text,
)
from . import data_path, load_parameter_from_ffitp


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


load_parameter_from_ffitp("forcefield.itp", str(data_path("oplsaam")))
protein_text = data_path("oplsaam", "protein.mol2").read_text()
# The reference OPLS-AA/M mol2 ships one duplicated NPRO hydrogen label.
# Original Xponge tolerates that input; normalize it for the stricter C++ registry.
protein_text = protein_text.replace(
    "  894   H2    0.000    0.000    0.000 oplsm_290    54     NPRO   0.330000",
    "  894   H1    0.000    0.000    0.000 oplsm_290    54     NPRO   0.330000",
)
register_residue_templates_from_mol2_text(protein_text)
_configure_protein_terminals()
