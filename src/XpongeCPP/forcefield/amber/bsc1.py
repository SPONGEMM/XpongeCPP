"""Register Amber bsc1 DNA templates and parameters."""

from ... import (
    configure_residue_template_head,
    configure_residue_template_tail,
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_pdb_residue_name_mapping,
    register_residue_templates_from_mol2_file,
)
from . import data_path

register_amber_parmdat_file(str(data_path("parm10.dat")))
register_amber_frcmod_file(str(data_path("parmbsc1.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("RNA.mol2")))
register_residue_templates_from_mol2_file(str(data_path("bsc1.mol2")))

for base in "ATCG":
    residue = f"D{base}"
    residue5 = f"{residue}5"
    residue3 = f"{residue}3"
    configure_residue_template_tail(residue, "O3'", 1.5, "C3'")
    configure_residue_template_tail(residue5, "O3'", 1.5, "C3'")
    configure_residue_template_head(residue, "P", 1.5, "OP2")
    configure_residue_template_head(residue3, "P", 1.5, "OP2")
    register_pdb_residue_name_mapping("head", residue, residue5)
    register_pdb_residue_name_mapping("tail", residue, residue3)
