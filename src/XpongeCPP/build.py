"""Legacy-style build module shim for XpongeCPP."""

from ._compat.workflows import build_bonded_force, get_mindsponge_system_energy

from ._compat.imports import reexport_module

reexport_module(
    "XpongeCPP",
    globals(),
    public=[
        "save_gro",
        "save_mol2",
        "save_pdb",
        "save_sponge_input",
        "save_sponge_input_raw",
        "Save_GRO",
        "Save_Mol2",
        "Save_PDB",
        "Save_SPONGE_Input",
        "Save_Sponge_Input",
    ],
)

__all__ = [
    "save_gro",
    "save_mol2",
    "save_pdb",
    "save_sponge_input",
    "save_sponge_input_raw",
    "Save_GRO",
    "Save_Mol2",
    "Save_PDB",
    "Save_SPONGE_Input",
    "Save_Sponge_Input",
    "build_bonded_force",
    "get_mindsponge_system_energy",
]
