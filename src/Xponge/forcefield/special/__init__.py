"""Legacy Xponge.forcefield.special shim."""

from XpongeCPP._compat.imports import extend_package_path

extend_package_path(globals(), "XpongeCPP.forcefield.special")

from . import fep, gb, min

__all__ = ["gb", "fep", "min"]
