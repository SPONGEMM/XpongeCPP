"""Compatibility shim around the shared PySCF QM backend."""

from __future__ import annotations

from ..qm.backends.pyscf_backend import require_numpy_pyscf
from ..qm.scheduler import compute_esp_on_grid as _compute_esp_on_grid
from ..qm.scheduler import run_scf as _run_scf


def build_backend_payload(assign, basis, charge, spin, opt, return_timings=False):
    scf_result = _run_scf(
        assign,
        backend="pyscf",
        basis=basis,
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


def compute_esp_on_grid(payload, grid_points_bohr):
    esp_result = _compute_esp_on_grid(payload["scf_result"], grid_points_bohr)
    return esp_result.electronic_esp_au
