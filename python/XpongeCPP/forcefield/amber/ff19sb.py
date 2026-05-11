"""Register Amber ff19SB templates, parameters, and CMAP terms."""

from ... import (
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
)
from . import data_path

register_amber_parmdat_file(str(data_path("parm19.dat")))
register_amber_frcmod_file(str(data_path("ff19SB.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("ff19SB.mol2")))
