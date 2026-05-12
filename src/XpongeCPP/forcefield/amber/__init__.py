"""Amber force-field compatibility namespace."""

from .. import package_data_path


def data_path(filename):
    return package_data_path("amber", filename)
