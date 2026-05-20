"""Backend-neutral RESP numerical helpers.

This module contains the current pure-Python RESP grid generation and fitting
logic. It is intentionally separated from any QM backend so the algorithm can
be tested independently and later migrated to C++.
"""

from __future__ import annotations

import math
import time

from .._core import fit_resp_from_esp_cpp as _fit_resp_from_esp_cpp
from .._core import fit_resp_from_esp_cpp_debug as _fit_resp_from_esp_cpp_debug
from .._core import generate_resp_mk_grid as _generate_resp_mk_grid_cpp


default_radius = {
    "H": 1.2,
    "C": 1.5,
    "N": 1.5,
    "O": 1.4,
    "P": 1.8,
    "S": 1.75,
    "F": 1.35,
    "Cl": 1.7,
    "Br": 2.3,
}


def _require_numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError("NumPy is required for RESP charge calculation") from exc
    return np


def fibonacci_grid(npoints, center, radius):
    if npoints <= 0:
        return []
    out = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(npoints):
        y = 1.0 - (2.0 * i + 1.0) / npoints
        r = math.sqrt(max(1.0 - y * y, 0.0))
        theta = golden_angle * i
        out.append(
            [
                center[0] + radius * math.cos(theta) * r,
                center[1] + radius * y,
                center[2] + radius * math.sin(theta) * r,
            ]
        )
    return out


def get_mk_grid(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
    np = _require_numpy()
    grids = []
    factor = area_density * 0.52918 * 0.52918 * 4 * np.pi
    real_radius = dict(default_radius)
    if radius:
        real_radius.update(radius)
    radial_layers = np.array([1.4 + 0.2 * i for i in range(layer)])
    for i, atom in enumerate(assign.atoms):
        if atom not in real_radius:
            raise KeyError(f"Radius for element {atom} not found")
        r0 = real_radius[atom] / 0.52918
        for r in r0 * radial_layers:
            grids.extend(fibonacci_grid(int(factor * r * r), atom_coordinates_bohr[i], r))
    grids = np.array(grids).reshape(-1, 3)
    for i, atom in enumerate(assign.atoms):
        r0 = 1.39 * real_radius[atom] / 0.52918
        t = np.linalg.norm(grids - atom_coordinates_bohr[i], axis=1)
        grids = grids[t >= r0, :]
    return grids


def get_mk_grid_cpp(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
    if radius is None:
        radius = {}
    return _generate_resp_mk_grid_cpp(list(assign.atoms), atom_coordinates_bohr, area_density, layer, radius)


def force_equivalence_q(q, extra_equivalence):
    np = _require_numpy()
    for eq_group in extra_equivalence:
        q_mean = np.mean(q[eq_group])
        q[eq_group] = q_mean
    return [float(x) for x in np.asarray(q).reshape(-1)]


def atom_judge(assign, atom, mask):
    element = "".join(ch for ch in mask if not ch.isdigit())
    digits = "".join(ch for ch in mask if ch.isdigit())
    if digits:
        return assign.atoms[atom] == element and len(assign.bonds[atom]) == int(digits)
    return assign.atoms[atom] == element


def find_tofit_second(assign, atom_count):
    tofit_second = []
    fit_group = {i: -1 for i in range(atom_count)}
    sublength = 0
    for i in range(atom_count):
        if atom_judge(assign, i, "C4"):
            fit_group[i] = len(tofit_second)
            tofit_second.append([i])
            hydrogens = [j for j in assign.bonds[i] if assign.atoms[j] == "H"]
            if hydrogens:
                for j in hydrogens:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(hydrogens)
                sublength += len(hydrogens) - 1
        if atom_judge(assign, i, "C3"):
            hydrogens = [j for j in assign.bonds[i] if assign.atoms[j] == "H"]
            if len(hydrogens) == 2:
                fit_group[i] = len(tofit_second)
                tofit_second.append([i])
                for j in hydrogens:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(hydrogens)
                sublength += 1
    return tofit_second, fit_group, sublength


def correct_extra_equivalence(tofit_second, fit_group, sublength, extra_equivalence, atom_numbers):
    if not extra_equivalence:
        return tofit_second, fit_group, sublength
    equi_group = []
    for eq in extra_equivalence:
        group = sorted({fit_group[eq_atom] for eq_atom in eq if fit_group[eq_atom] != -1})
        equi_group.append(group)
    all_groups_list = sorted({fit_group[atom] for atom in range(atom_numbers)})
    group_map = {i: i for i in all_groups_list}
    for eq in equi_group:
        for group in eq:
            group_map[group] = eq[0]
    temp_max = 0
    for group in all_groups_list:
        if group == -1:
            continue
        if group_map[group] == group:
            group_map[group] = temp_max
            temp_max += 1
        else:
            group_map[group] = group_map[group_map[group]]
    old = tofit_second
    tofit_second = [[] for _ in range(temp_max)]
    for i, group in enumerate(old):
        tofit_second[group_map[i]].extend(group)
        sublength -= len(group) - 1
    for group in tofit_second:
        sublength += len(group) - 1
    for atom in range(atom_numbers):
        fit_group[atom] = group_map[fit_group[atom]]
    return tofit_second, fit_group, sublength


def resp_scf_kernel(assign, atom_count, a, b, matrix_a, matrix_a0, matrix_b, q):
    np = _require_numpy()
    q = np.asarray(q, dtype=float).reshape(-1)
    step = 0
    q_last_step = q
    while step == 0 or np.max(np.abs(q - q_last_step)) > 1e-4:
        step += 1
        q_last_step = q
        for i in range(atom_count):
            if assign.atoms[i] != "H":
                matrix_a[i][i] = matrix_a0[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
        q = np.dot(np.linalg.inv(matrix_a), matrix_b).reshape(-1)[:-1]
    return q


def get_a20_and_b20(total_length, tofit_second, fit_group, sublength, atom_count, matrix_a0, matrix_b, charge, q):
    np = _require_numpy()
    a20 = np.zeros((total_length, total_length))
    count = len(tofit_second)
    for i in range(atom_count):
        if fit_group[i] == -1:
            fit_group[i] = count
            count += 1
        a20[atom_count - sublength][fit_group[i]] += 1
        a20[fit_group[i]][atom_count - sublength] += 1
    b20 = np.zeros(total_length)
    for i in range(atom_count):
        b20[fit_group[i]] += matrix_b[i]
        for j in range(atom_count):
            a20[fit_group[i]][fit_group[j]] += matrix_a0[i][j]
    b20[atom_count - sublength] = charge
    count = 0
    for i in range(atom_count):
        if fit_group[i] >= len(tofit_second):
            b20[atom_count - sublength + count + 1] = q[i]
            a20[atom_count - sublength + count + 1][len(tofit_second) + count] = 1
            a20[len(tofit_second) + count][atom_count - sublength + count + 1] = 1
            count += 1
    return a20, b20


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
    )["final_charges"]


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
    np = _require_numpy()
    if extra_equivalence is None:
        extra_equivalence = []
    atom_count = len(assign.atoms)
    timings = {}
    total_start = time.perf_counter()
    assembly_start = time.perf_counter()
    vnuc = 0
    matrix_a0 = np.zeros((atom_count, atom_count))
    for i in range(atom_count):
        r = atom_coordinates_bohr[i]
        z = nuclear_charges[i]
        rp = r - grid_points_bohr
        for j in range(atom_count):
            rpj = atom_coordinates_bohr[j] - grid_points_bohr
            matrix_a0[i][j] = np.sum(1.0 / np.linalg.norm(rp, axis=1) / np.linalg.norm(rpj, axis=1))
        vnuc += z / np.einsum("xi,xi->x", rp, rp) ** 0.5
    matrix_a0 = np.hstack((matrix_a0, np.ones(atom_count).reshape(-1, 1)))
    temp = np.ones(atom_count + 1)
    temp[-1] = 0
    matrix_a0 = np.vstack((matrix_a0, temp.reshape(1, -1)))
    matrix_a = np.zeros_like(matrix_a0)
    matrix_a[:] = matrix_a0

    mep = vnuc - esp_values_au
    matrix_b = np.zeros((atom_count + 1))
    for i in range(atom_count):
        r = atom_coordinates_bohr[i]
        rp = np.linalg.norm(r - grid_points_bohr, axis=1)
        matrix_b[i] = np.sum(mep / rp)
    matrix_b[-1] = charge
    matrix_b = matrix_b.reshape(-1, 1)
    matrix_b_flat = matrix_b.reshape(-1)
    timings["assembly"] = time.perf_counter() - assembly_start

    stage1_start = time.perf_counter()
    q = np.dot(np.linalg.inv(matrix_a), matrix_b).reshape(-1)[:-1]
    esp_charges = force_equivalence_q(q.copy(), extra_equivalence)

    if only_esp:
        timings["stage1"] = time.perf_counter() - stage1_start
        timings["stage2"] = 0.0
        timings["total"] = time.perf_counter() - total_start
        return {
            "esp_charges": esp_charges,
            "stage1_charges": esp_charges,
            "final_charges": esp_charges,
            "timings": timings,
        }

    q = resp_scf_kernel(assign, atom_count, a1, 0.1, matrix_a, matrix_a0, matrix_b_flat, q)
    stage1_charges = force_equivalence_q(q.copy(), extra_equivalence)
    timings["stage1"] = time.perf_counter() - stage1_start
    if not two_stage:
        timings["stage2"] = 0.0
        timings["total"] = time.perf_counter() - total_start
        return {
            "esp_charges": esp_charges,
            "stage1_charges": stage1_charges,
            "final_charges": stage1_charges,
            "timings": timings,
        }

    stage2_start = time.perf_counter()
    tofit_second, fit_group, sublength = find_tofit_second(assign, atom_count)
    tofit_second, fit_group, sublength = correct_extra_equivalence(
        tofit_second, fit_group, sublength, extra_equivalence, atom_count
    )
    if tofit_second:
        total_length = atom_count - sublength + 1 + atom_count - sublength - len(tofit_second)
        a20, b20 = get_a20_and_b20(
            total_length, tofit_second, fit_group, sublength, atom_count, matrix_a0, matrix_b_flat, charge, q
        )
        matrix_a = np.zeros_like(a20)
        matrix_a[:] = a20[:]
        matrix_b = b20.reshape(-1, 1)
        q_temp = np.dot(np.linalg.inv(matrix_a), matrix_b).reshape(-1)[:-1]
        step = 0
        q_last_step = q_temp
        while step == 0 or np.max(np.abs(q_temp - q_last_step)) > 1e-4:
            step += 1
            q_last_step = q_temp
            for i in range(atom_count - sublength):
                if assign.atoms[i] != "H":
                    matrix_a[i][i] = a20[i][i] + a2 / np.sqrt(q_last_step[i] * q_last_step[i] + 0.1 * 0.1)
            q_temp = np.dot(np.linalg.inv(matrix_a), matrix_b).reshape(-1)[:-1]
        for i, group in enumerate(tofit_second):
            for j in group:
                q[j] = q_temp[i]
    final_charges = force_equivalence_q(q, extra_equivalence)
    timings["stage2"] = time.perf_counter() - stage2_start
    timings["total"] = time.perf_counter() - total_start
    return {
        "esp_charges": esp_charges,
        "stage1_charges": stage1_charges,
        "final_charges": final_charges,
        "timings": timings,
    }


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
