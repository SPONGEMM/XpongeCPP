"""Legacy charge_base compatibility helpers."""

from __future__ import annotations

from ...helper import AtomType, Molecule, set_global_alternative_names

AtomType.Add_Property({"charge": float})


@Molecule.Set_MindSponge_Todo("charge")
def _do_charge(self, sys_kwarg, ene_kwarg, use_pbc):
    sys_kwarg.setdefault("atom_charge", []).append([float(getattr(atom, "charge", 0.0)) for atom in self.atoms])

    if "charge" in ene_kwarg:
        return

    def _build_coulomb_energy(system, ene_kwarg):
        del ene_kwarg
        from mindsponge.potential import CoulombEnergy

        kwargs = {
            "atom_charge": system.atom_charge,
            "length_unit": "A",
            "energy_unit": "kcal/mol",
        }
        if hasattr(system, "pbc_box") and use_pbc:
            kwargs["pbc_box"] = system.pbc_box
        return CoulombEnergy(**kwargs)

    ene_kwarg["charge"] = {"function": _build_coulomb_energy}


AtomType.Set_Property_Unit("charge", "charge", "e")

set_global_alternative_names()

__all__ = []
