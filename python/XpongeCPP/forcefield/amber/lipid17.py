"""Register Amber Lipid17 templates and parameters."""

from ... import register_amber_parmdat_file, register_residue_templates_from_mol2_file
from . import data_path

register_amber_parmdat_file(str(data_path("lipid17.dat")))
register_residue_templates_from_mol2_file(str(data_path("lipid17.mol2")))
