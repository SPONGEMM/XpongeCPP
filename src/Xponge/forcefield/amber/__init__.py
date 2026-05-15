"""Legacy Xponge.forcefield.amber shim."""

from XpongeCPP._compat.imports import extend_package_path, reexport_module

reexport_module(
    "XpongeCPP.forcefield.amber",
    globals(),
    public=["data_path", "load_parameters_from_frcmod", "Load_Parameters_From_Frcmod"],
)
extend_package_path(globals(), "XpongeCPP.forcefield.amber")
