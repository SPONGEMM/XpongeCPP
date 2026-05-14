"""Register SPC/E water and Amber ion templates from packaged data."""

from ... import (
    register_amber_bond_parameter,
    register_amber_frcmod_file,
    register_amber_lj_parameter,
    register_template_molecule_from_mol2_file,
)
from . import data_path

register_template_molecule_from_mol2_file(str(data_path("spce.mol2")))
register_amber_bond_parameter("OW", "HW", 553.0, 1.0000)
register_amber_bond_parameter("HW", "HW", 553.0, 1.6330)
register_amber_lj_parameter("OW", "OW", 0.1553, 1.7767)
register_amber_lj_parameter("HW", "HW", 0.0, 0.0)

register_amber_frcmod_file(str(data_path("ions1lm_126_spce.frcmod")))
register_amber_frcmod_file(str(data_path("ionsjc_spce.frcmod")))
register_amber_frcmod_file(str(data_path("ions234lm_126_spce.frcmod")))
register_template_molecule_from_mol2_file(str(data_path("atomic_ions.mol2")))
