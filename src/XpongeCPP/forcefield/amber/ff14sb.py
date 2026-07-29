"""Register ff14SB templates from packaged Amber force-field data."""

from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("protein", "ff14sb")

from ... import (
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_ff14sb,
    register_residue_templates_from_mol2_file,
)
from . import configure_proline_like_terminal_mapping, data_path

register_ff14sb()
register_amber_parmdat_file(str(data_path("parm10.dat")))
register_amber_frcmod_file(str(data_path("ff14SB.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("ff14SB.mol2")))
configure_proline_like_terminal_mapping("HYP", "CHYP")
