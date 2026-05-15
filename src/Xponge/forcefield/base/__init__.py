"""Legacy Xponge.forcefield.base shim."""

from XpongeCPP._compat.imports import extend_package_path

extend_package_path(globals(), "XpongeCPP.forcefield.base")

from . import angle_base, bond_base, charge_base, cmap_base, dihedral_base, exclude_base, lj_base, mass_base
from . import nb14_base, virtual_atom_base

__all__ = [
    "angle_base",
    "bond_base",
    "charge_base",
    "cmap_base",
    "dihedral_base",
    "exclude_base",
    "lj_base",
    "mass_base",
    "nb14_base",
    "virtual_atom_base",
]
