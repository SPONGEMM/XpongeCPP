"""Register ff14SB templates from packaged Amber force-field data."""

from ... import (
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_ff14sb,
    register_residue_templates_from_mol2_file,
)
from . import data_path

register_ff14sb()
register_amber_parmdat_file(str(data_path("parm10.dat")))
register_amber_frcmod_file(str(data_path("ff14SB.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("ff14SB.mol2")))
