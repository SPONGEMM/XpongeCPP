"""Register ff14SB templates from packaged Amber force-field data."""

from ... import register_ff14sb, register_residue_templates_from_mol2_file
from . import data_path

register_ff14sb()
register_residue_templates_from_mol2_file(str(data_path("ff14SB.mol2")))
