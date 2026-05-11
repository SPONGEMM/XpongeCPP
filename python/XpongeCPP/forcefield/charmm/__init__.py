"""CHARMM force-field compatibility namespace."""

from importlib import resources


def data_path(*parts):
    return resources.files("XpongeCPP").joinpath("data", "reference_forcefield", "charmm", *parts)
