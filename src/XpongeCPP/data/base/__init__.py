"""Compatibility bridge for packaged legacy data modules.

Legacy data modules under ``XpongeCPP.data.amber`` still import ``..base`` using
the original Xponge package layout.  Re-export the forcefield base shims here so
those relative imports keep working.
"""

from ...forcefield.base import angle_base, bond_base, charge_base, exclude_base, lj_base, mass_base
from ...forcefield.base import (
    cmap_base,
    dihedral_base,
    nb14_base,
    virtual_atom_base,
)

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
