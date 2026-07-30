"""Legacy Xponge.forcefield.edip shim."""

from XpongeCPP._compat.imports import extend_package_path, reexport_module

reexport_module("XpongeCPP.forcefield.edip", globals(), public=["data_path"])
extend_package_path(globals(), "XpongeCPP.forcefield.edip")
