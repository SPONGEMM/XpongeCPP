"""Register Amber OL3 RNA templates."""

from ... import register_amber_parmdat_file, register_residue_templates_from_mol2_file
from . import data_path

register_amber_parmdat_file(str(data_path("parm10.dat")))
register_residue_templates_from_mol2_file(str(data_path("RNA.mol2")))
