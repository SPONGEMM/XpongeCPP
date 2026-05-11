"""Register packaged CHARMM36 GROMACS-port parameters."""

from ... import load_gromacs_topology_file
from . import data_path

load_gromacs_topology_file(str(data_path("charmm36", "forcefield.itp")))
