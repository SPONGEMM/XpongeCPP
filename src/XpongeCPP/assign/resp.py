"""RESP orchestration module.

The long-term target is a Psi4 backend plus a C++ RESP core. During migration,
this module stays as the stable Python compatibility entrypoint and delegates
the current behavior to a backend module plus a backend-neutral Python core.
"""

from __future__ import annotations

import sys

from ..qm import scheduler as qm_scheduler
from ..qm.resp_basis import resolve_default_resp_basis, resolve_resp_basis
from ..qm.resp_parameters import get_resp_radius_overrides, normalize_element_symbol
from . import psi4_backend, pyscf_backend, resp_core


_BACKEND_MODULES = {
    "pyscf": pyscf_backend,
    "psi4": psi4_backend,
}
_CORE_MODULES = {"python", "cpp"}
_RESP_BASE_REFERENCES = (
    "Bayly1993_RESP",
    "SinghKollman1984_MK",
    "BeslerMerzKollman1990_ESP",
)


def _resolve_basis_for_resp(assign, basis):
    elements = set(assign.atoms)
    if basis is None:
        return resolve_default_resp_basis(elements)
    return resolve_resp_basis(basis, elements)


def _merge_resp_radii(assign, radius):
    resolved = get_resp_radius_overrides(assign.atoms)
    if radius:
        resolved.update({normalize_element_symbol(element): float(value) for element, value in dict(radius).items()})
    return resolved


def _build_resp_metadata(assign, resolved_basis, radius):
    references = []
    for reference in _RESP_BASE_REFERENCES + tuple(resolved_basis.references):
        if reference not in references:
            references.append(reference)
    return {
        "method": "RESP",
        "esp_setup": "standard_mk_hf",
        "basis_family": resolved_basis.label,
        "basis": resolved_basis.basis,
        "ecp": resolved_basis.ecp,
        "cart": resolved_basis.cart,
        "radius_set": "standard MK radii",
        "radius": dict(radius),
        "references": references,
    }


def get_resp_setup_metadata(assign, basis=None, radius=None):
    resolved_basis = _resolve_basis_for_resp(assign, basis)
    resolved_radius = _merge_resp_radii(assign, radius)
    return _build_resp_metadata(assign, resolved_basis, resolved_radius)


def _build_backend_payload(backend_module, assign, resolved_basis, charge, spin, opt):
    try:
        return backend_module.build_backend_payload(
            assign,
            resolved_basis.basis,
            charge,
            spin,
            opt,
            ecp=resolved_basis.ecp,
            cart=resolved_basis.cart,
        )
    except TypeError as exc:
        if "unexpected keyword" not in str(exc):
            raise
        return backend_module.build_backend_payload(assign, resolved_basis.basis, charge, spin, opt)


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


def resp_fit(assign, basis=None, opt=False, charge=None, spin=0, extra_equivalence=None,
             grid_density=6, grid_cell_layer=4, radius=None, a1=0.0005, a2=0.001,
             two_stage=True, only_esp=False, backend=None, core=None,
             esp_memory_limit=None, esp_chunk_policy="auto", esp_safety_factor=0.8,
             return_metadata=False):
    if extra_equivalence is None:
        extra_equivalence = []
    if charge is None:
        charge = int(round(sum(assign.charges)))
    backend_name = _normalize_backend_name(backend)
    _normalize_core_name(core)
    backend_module = _BACKEND_MODULES[backend_name]
    resolved_basis = _resolve_basis_for_resp(assign, basis)
    resolved_radius = _merge_resp_radii(assign, radius)
    metadata = _build_resp_metadata(assign, resolved_basis, resolved_radius)
    try:
        payload = _build_backend_payload(backend_module, assign, resolved_basis, charge, spin, opt)
    except ImportError as exc:
        _backend_import_or_hint(backend_name, exc)
    grids = resp_core.get_mk_grid(
        assign,
        payload["atom_coordinates_bohr"],
        area_density=grid_density,
        layer=grid_cell_layer,
        radius=resolved_radius,
    )
    electron_esp = backend_module.compute_esp_on_grid(
        payload,
        grids,
        memory_limit=esp_memory_limit,
        chunk_policy=esp_chunk_policy,
        safety_factor=esp_safety_factor,
    )
    charges = resp_core.fit_resp_from_esp(
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
    if return_metadata:
        return {"charges": charges, "metadata": metadata}
    return charges
