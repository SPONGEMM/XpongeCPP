"""Register Amber GAFF parameters from packaged data."""

from ... import register_amber_parmdat_file
from . import data_path

register_amber_parmdat_file(str(data_path("gaff.dat")))
