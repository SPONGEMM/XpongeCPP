"""Register packaged OPLS-AA/M GROMACS parameters."""

from ... import load_opls_itp_file
from . import data_path

load_opls_itp_file(str(data_path("oplsaam", "forcefield.itp")))
