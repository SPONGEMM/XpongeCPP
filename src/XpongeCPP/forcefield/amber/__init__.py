"""Amber force-field compatibility namespace."""

from .. import package_data_path
from ... import register_amber_frcmod_file


def data_path(filename):
    return package_data_path("amber", filename)


def load_parameters_from_frcmod(filename, prefix=True):
    return register_amber_frcmod_file(str(filename))


Load_Parameters_From_Frcmod = load_parameters_from_frcmod
