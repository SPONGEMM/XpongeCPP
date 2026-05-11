"""Register Amber OL15 DNA templates and parameter updates."""

from ... import (
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
)
from . import data_path

register_amber_parmdat_file(str(data_path("parm10.dat")))
register_amber_frcmod_file(str(data_path("OL15.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("ol15.mol2")))
