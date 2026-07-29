"""Legacy-style load module shim for XpongeCPP."""

from ._compat.imports import reexport_module

reexport_module(
    "XpongeCPP",
    globals(),
    public=[
        "GromacsTopologyIterator",
        "load_coordinate",
        "load_ffitp",
        "load_frcmod",
        "load_gro",
        "load_mol2",
        "load_molitp",
        "load_molpsf",
        "load_mmcif",
        "load_pdb",
        "load_rst7",
    ],
)
