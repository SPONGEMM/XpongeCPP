"""Register TIP3P water and Amber ion templates from packaged data."""

from ... import (
    AtomType,
    register_amber_frcmod_file,
    register_residue_template_alias,
    register_residue_templates_from_mol2_file,
    register_tip3p,
)
from . import data_path

AtomType.New_From_String(
    """
name mass    charge[e]  LJtype
HW   1.008    0.417       HW
OW   16      -0.834       OW
"""
)

register_tip3p()
register_residue_templates_from_mol2_file(str(data_path("tip3p.mol2")))
register_residue_template_alias("H2O", "WAT")
register_residue_template_alias("HOH", "WAT")
register_amber_frcmod_file(str(data_path("ions1lm_126_tip3p.frcmod")))
register_amber_frcmod_file(str(data_path("ionsjc_tip3p.frcmod")))
register_amber_frcmod_file(str(data_path("ions234lm_126_tip3p.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("atomic_ions.mol2")))
