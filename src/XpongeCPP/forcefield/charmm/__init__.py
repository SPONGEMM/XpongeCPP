"""CHARMM force-field compatibility namespace."""

from .. import repository_reference_forcefield_path
from ... import set_lj_combining_rule


set_lj_combining_rule("lorentz_berthelot")


def data_path(*parts):
    return repository_reference_forcefield_path("charmm", *parts)
