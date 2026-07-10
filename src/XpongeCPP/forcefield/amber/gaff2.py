"""Register Amber GAFF2 parameters from packaged data."""

from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("small_molecule", "gaff2")

from ... import AtomType, implemented_gaff2_assign_types, register_amber_parmdat_file
from . import data_path

register_amber_parmdat_file(str(data_path("gaff2.dat")))
AtomType.New_From_String("\n".join(implemented_gaff2_assign_types()))
