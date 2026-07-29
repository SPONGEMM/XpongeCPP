"""Compatibility shim around the shared Psi4 QM backend."""

from __future__ import annotations

from ..qm.backends.psi4_backend import ANGSTROM_PER_BOHR, require_numpy_psi4
from ..qm.scheduler import compute_esp_on_grid as _compute_esp_on_grid
from ..qm.scheduler import run_scf as _run_scf


def build_backend_payload(assign, basis, charge, spin, opt, return_timings=False, ecp=None, cart=None):
    scf_result = _run_scf(
        assign,
        backend="psi4",
        basis=basis,
        ecp=ecp,
        cart=cart,
        charge=charge,
        spin=spin,
        optimize_geometry=opt,
        return_timings=return_timings,
    )
    payload = {
        "atom_symbols": scf_result.atom_symbols,
        "atom_coordinates_bohr": scf_result.coordinates_bohr,
        "nuclear_charges": scf_result.nuclear_charges,
        "charge": scf_result.charge,
        "spin": scf_result.spin,
        "scf_result": scf_result,
    }
    if return_timings:
        payload["timings"] = dict(scf_result.timings)
    return payload


def compute_esp_on_grid(payload, grid_points_bohr, *, memory_limit=None, chunk_policy="auto", safety_factor=0.8):
    esp_result = _compute_esp_on_grid(
        payload["scf_result"],
        grid_points_bohr,
        memory_limit=memory_limit,
        chunk_policy=chunk_policy,
        safety_factor=safety_factor,
    )
    return esp_result.electronic_esp_au
