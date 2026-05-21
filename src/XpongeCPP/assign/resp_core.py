"""RESP numerical wrappers.

Python now keeps only a thin orchestration layer for RESP. The numerical grid
generation and charge fitting paths are delegated to the C++ core.
"""

from __future__ import annotations

from .._core import fit_resp_from_esp_cpp as _fit_resp_from_esp_cpp
from .._core import fit_resp_from_esp_cpp_debug as _fit_resp_from_esp_cpp_debug
from .._core import generate_resp_mk_grid as _generate_resp_mk_grid_cpp


def get_mk_grid(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
    """Compatibility wrapper that routes MK-grid generation to C++."""
    if radius is None:
        radius = {}
    return _generate_resp_mk_grid_cpp(list(assign.atoms), atom_coordinates_bohr, area_density, layer, radius)


def get_mk_grid_cpp(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
    """Explicit C++ entry kept for compatibility with older tests and scripts."""
    return get_mk_grid(
        assign,
        atom_coordinates_bohr,
        area_density=area_density,
        layer=layer,
        radius=radius,
    )


def fit_resp_from_esp(
    assign,
    atom_coordinates_bohr,
    nuclear_charges,
    grid_points_bohr,
    esp_values_au,
    charge,
    extra_equivalence=None,
    a1=0.0005,
    a2=0.001,
    two_stage=True,
    only_esp=False,
):
    """Compatibility wrapper that routes RESP fitting to C++."""
    if extra_equivalence is None:
        extra_equivalence = []
    return _fit_resp_from_esp_cpp(
        assign,
        atom_coordinates_bohr,
        nuclear_charges,
        grid_points_bohr,
        esp_values_au,
        int(charge),
        extra_equivalence,
        a1,
        a2,
        two_stage,
        only_esp,
    )


def fit_resp_from_esp_debug(
    assign,
    atom_coordinates_bohr,
    nuclear_charges,
    grid_points_bohr,
    esp_values_au,
    charge,
    extra_equivalence=None,
    a1=0.0005,
    a2=0.001,
    two_stage=True,
    only_esp=False,
):
    """Compatibility wrapper that routes RESP debug fitting to C++."""
    if extra_equivalence is None:
        extra_equivalence = []
    return _fit_resp_from_esp_cpp_debug(
        assign,
        atom_coordinates_bohr,
        nuclear_charges,
        grid_points_bohr,
        esp_values_au,
        int(charge),
        extra_equivalence,
        a1,
        a2,
        two_stage,
        only_esp,
    )


def fit_resp_from_esp_cpp(
    assign,
    atom_coordinates_bohr,
    nuclear_charges,
    grid_points_bohr,
    esp_values_au,
    charge,
    extra_equivalence=None,
    a1=0.0005,
    a2=0.001,
    two_stage=True,
    only_esp=False,
):
    return fit_resp_from_esp(
        assign,
        atom_coordinates_bohr,
        nuclear_charges,
        grid_points_bohr,
        esp_values_au,
        charge,
        extra_equivalence=extra_equivalence,
        a1=a1,
        a2=a2,
        two_stage=two_stage,
        only_esp=only_esp,
    )


def fit_resp_from_esp_cpp_debug(
    assign,
    atom_coordinates_bohr,
    nuclear_charges,
    grid_points_bohr,
    esp_values_au,
    charge,
    extra_equivalence=None,
    a1=0.0005,
    a2=0.001,
    two_stage=True,
    only_esp=False,
):
    return fit_resp_from_esp_debug(
        assign,
        atom_coordinates_bohr,
        nuclear_charges,
        grid_points_bohr,
        esp_values_au,
        charge,
        extra_equivalence=extra_equivalence,
        a1=a1,
        a2=a2,
        two_stage=two_stage,
        only_esp=only_esp,
    )
