"""Register Amber ff19SB templates, parameters, and CMAP terms."""

from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("protein", "ff19sb")

from ... import (
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
)
from . import configure_proline_like_terminal_mapping, data_path

register_amber_parmdat_file(str(data_path("parm19.dat")))
register_amber_frcmod_file(str(data_path("ff19SB.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("ff19SB.mol2")))
register_residue_templates_from_mol2_file(str(data_path("ff19SB_nhyp.mol2")))
configure_proline_like_terminal_mapping("HYP", "CHYP", "NHYP")
