"""Register TIP4P water templates from packaged data."""

from ... import (
    register_amber_bond_parameter,
    register_amber_lj_parameter,
    register_template_molecule_from_mol2_file,
    register_template_virtual_atom2,
)
from . import data_path

register_template_molecule_from_mol2_file(str(data_path("tip4pew.mol2")))
register_amber_bond_parameter("OW", "HW", 553.0, 0.9572)
register_amber_bond_parameter("HW", "HW", 553.0, 1.5136)
register_amber_lj_parameter("OW", "OW", 0.155, 3.154 / 2.0)
register_amber_lj_parameter("HW", "HW", 0.0, 0.0)
register_amber_lj_parameter("EP", "EPW", 0.0, 0.0)
register_template_virtual_atom2("WAT", "EPW", "O", "H1", "H2", 0.1279696, 0.1279696)
