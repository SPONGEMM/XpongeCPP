"""EDIP force-field compatibility namespace."""

from .. import repository_reference_forcefield_path


def data_path(*parts):
    return repository_reference_forcefield_path("edip", *parts)
