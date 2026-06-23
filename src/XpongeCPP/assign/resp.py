"""RESP orchestration module.

The long-term target is a Psi4 backend plus a C++ RESP core. During migration,
this module stays as the stable Python compatibility entrypoint and delegates
the current behavior to a backend module plus a backend-neutral Python core.
"""

from __future__ import annotations

import sys

from ..qm import scheduler as qm_scheduler
from . import psi4_backend, pyscf_backend, resp_core


_BACKEND_MODULES = {
    "pyscf": pyscf_backend,
    "psi4": psi4_backend,
}
_CORE_MODULES = {"python", "cpp"}


def _normalize_backend_name(backend):
    try:
        return qm_scheduler.normalize_backend_name(backend)
    except ValueError as exc:
        raise ValueError(str(exc).replace("QM backend", "RESP backend", 1)) from exc


def _backend_import_or_hint(backend_name, exc):
    message = str(exc)
    if backend_name == "pyscf" and sys.platform.startswith("win"):
        message += " On Windows, install Psi4 via conda-forge or the official Psi4 installer and call calculate_charge('resp', backend='psi4', ...)."
    elif backend_name == "psi4":
        message += " On Windows, Psi4 is not installed through pip by default; install it via conda-forge or the official Psi4 installer and retry."
    raise ImportError(message) from exc


def _normalize_core_name(core):
    if core is None:
        return "cpp"
    core_name = str(core).strip().lower()
    if core_name not in _CORE_MODULES:
        supported = ", ".join(sorted(_CORE_MODULES))
        raise ValueError(f"RESP core should be one of: {supported}")
    return "cpp"


def resp_fit(assign, basis="6-31g*", opt=False, charge=None, spin=0, extra_equivalence=None,
             grid_density=6, grid_cell_layer=4, radius=None, a1=0.0005, a2=0.001,
             two_stage=True, only_esp=False, backend=None, core=None,
             esp_memory_limit=None, esp_chunk_policy="auto", esp_safety_factor=0.8):
    if extra_equivalence is None:
        extra_equivalence = []
    if charge is None:
        charge = int(round(sum(assign.charges)))
    backend_name = _normalize_backend_name(backend)
    _normalize_core_name(core)
    backend_module = _BACKEND_MODULES[backend_name]
    try:
        payload = backend_module.build_backend_payload(assign, basis, charge, spin, opt)
    except ImportError as exc:
        _backend_import_or_hint(backend_name, exc)
    grids = resp_core.get_mk_grid(
        assign,
        payload["atom_coordinates_bohr"],
        area_density=grid_density,
        layer=grid_cell_layer,
        radius=radius,
    )
    electron_esp = backend_module.compute_esp_on_grid(
        payload,
        grids,
        memory_limit=esp_memory_limit,
        chunk_policy=esp_chunk_policy,
        safety_factor=esp_safety_factor,
    )
    return resp_core.fit_resp_from_esp(
        assign,
        atom_coordinates_bohr=payload["atom_coordinates_bohr"],
        nuclear_charges=payload["nuclear_charges"],
        grid_points_bohr=grids,
        esp_values_au=electron_esp,
        charge=charge,
        extra_equivalence=extra_equivalence,
        a1=a1,
        a2=a2,
        two_stage=two_stage,
        only_esp=only_esp,
    )
