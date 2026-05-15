"""Legacy mass_base compatibility helpers."""

from __future__ import annotations

from ...helper import AtomType, Molecule, set_global_alternative_names

AtomType.Add_Property({"mass": float})


@Molecule.Set_MindSponge_Todo("mass")
def _do_mass(self, sys_kwarg, ene_kwarg, use_pbc):
    del ene_kwarg
    del use_pbc
    sys_kwarg.setdefault("atom_mass", []).append([float(getattr(atom, "mass", 0.0)) for atom in self.atoms])


AtomType.Set_Property_Unit("mass", "mass", "amu")

set_global_alternative_names()

__all__ = []
