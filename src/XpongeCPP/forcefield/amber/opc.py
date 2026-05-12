"""Register OPC water and Amber ion templates from packaged data."""

from ... import (
    register_amber_bond_parameter,
    register_amber_frcmod_file,
    register_amber_lj_parameter,
    register_template_molecule_from_mol2_file,
    register_template_virtual_atom2,
)
from . import data_path

register_template_molecule_from_mol2_file(str(data_path("opc.mol2")))
register_amber_bond_parameter("OW", "HW", 553.0, 0.87243313)
register_amber_bond_parameter("HW", "HW", 553.0, 1.37120510)
register_amber_lj_parameter("OW", "OW", 0.2128008130, 1.777167268)
register_amber_lj_parameter("HW", "HW", 0.0, 0.0)
register_amber_lj_parameter("EP", "EPW", 0.0, 0.0)
register_template_virtual_atom2("WAT", "EPW", "O", "H1", "H2", 0.1477206, 0.1477206)

register_template_molecule_from_mol2_file(str(data_path("atomic_ions.mol2")))
register_amber_frcmod_file(str(data_path("ions1lm_126_tip4pew.frcmod")))
register_amber_frcmod_file(str(data_path("ionsjc_tip4pew.frcmod")))
register_amber_frcmod_file(str(data_path("ions234lm_126_tip4pew.frcmod")))
