"""Legacy Xponge.forcefield.opls shim."""

from XpongeCPP._compat.imports import extend_package_path, reexport_module

reexport_module("XpongeCPP.forcefield.opls", globals(), public=["data_path", "load_parameter_from_ffitp"])
extend_package_path(globals(), "XpongeCPP.forcefield.opls")
