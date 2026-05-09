"""Register TIP3P water and Amber ion templates from packaged data."""

from ... import register_residue_templates_from_mol2_file, register_tip3p
from . import data_path

register_tip3p()
register_residue_templates_from_mol2_file(str(data_path("tip3p.mol2")))
register_residue_templates_from_mol2_file(str(data_path("atomic_ions.mol2")))
