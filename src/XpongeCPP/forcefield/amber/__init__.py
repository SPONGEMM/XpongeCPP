"""Amber force-field compatibility namespace."""

from .. import package_data_path
from ... import register_amber_frcmod_file, register_amber_nb14_scale, set_lj_combining_rule


set_lj_combining_rule("lorentz_berthelot")


def data_path(filename):
    return package_data_path("amber", filename)


def load_parameters_from_frcmod(filename, prefix=True):
    set_lj_combining_rule("lorentz_berthelot")
    register_amber_nb14_scale("X", "X", 0.5, 0.833333)
    return register_amber_frcmod_file(str(filename))


Load_Parameters_From_Frcmod = load_parameters_from_frcmod
