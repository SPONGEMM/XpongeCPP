"""Amber force-field compatibility namespace."""

from importlib import resources


def data_path(filename):
    return resources.files("XpongeCPP").joinpath("data", "amber", filename)
