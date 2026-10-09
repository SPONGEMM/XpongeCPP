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

from ._nucleic import configure_nucleic_templates

configure_nucleic_templates(("DA", "DC", "DG", "DT"))
