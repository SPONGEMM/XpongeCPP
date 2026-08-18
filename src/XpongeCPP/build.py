"""Legacy-style build module shim for XpongeCPP."""

from .helper.math import guess_element_from_mass
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
        "save_sponge_input_bundle",
        "save_sponge_input_raw",
        "Save_GRO",
        "Save_Mol2",
        "Save_PDB",
        "Save_SPONGE_Input",
        "Save_Sponge_Input",
    ],
)


def _pdb_guess_element(atom):
    """Return the two-column PDB element token used by Xponge writers."""

    element = getattr(atom, "element", None)
    if not element:
        element = getattr(getattr(atom, "type", None), "element", None)
    if not element:
        mass = getattr(atom, "mass", None)
        if mass is not None and mass > 0:
            try:
                element = guess_element_from_mass(mass)
            except Exception:
                element = None
    if not element:
        name = getattr(atom, "name", "")
        letters = "".join(character for character in name if character.isalpha())
        if letters:
            element = (
                letters[0].upper() + letters[1].lower()
                if len(letters) >= 2 and letters[1].islower()
                else letters[0].upper()
            )
        else:
            element = "X"
    if len(element) == 1:
        return f"{element:>2}"
    return f"{element[:2]:>2}"


__all__ = [
    "save_gro",
    "save_mol2",
    "save_pdb",
    "save_sponge_input",
    "save_sponge_input_bundle",
    "save_sponge_input_raw",
    "Save_GRO",
    "Save_Mol2",
    "Save_PDB",
    "Save_SPONGE_Input",
    "Save_Sponge_Input",
    "build_bonded_force",
    "get_mindsponge_system_energy",
]
