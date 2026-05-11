"""Register GLYCAM-06j base carbohydrate templates and parameters."""

from .... import register_amber_parmdat_file, register_residue_templates_from_mol2_file
from .. import data_path

register_amber_parmdat_file(str(data_path("glycam_06j/GLYCAM_06j.dat")))
register_residue_templates_from_mol2_file(str(data_path("glycam_06j/terminal.mol2")))
