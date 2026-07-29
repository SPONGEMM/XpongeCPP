"""RESP numerical wrappers.

Python now keeps only a thin orchestration layer for RESP. The numerical grid
generation and charge fitting paths are delegated to the C++ core.
"""

from __future__ import annotations

import numpy as np

from .._core import fit_resp_from_esp_cpp as _fit_resp_from_esp_cpp
from .._core import fit_resp_from_esp_cpp_debug as _fit_resp_from_esp_cpp_debug
from .._core import generate_resp_mk_grid as _generate_resp_mk_grid_cpp
from .._core import (
    solve_resp_constrained_quadratic_cpp
    as _solve_resp_constrained_quadratic_cpp,
)


def _prepare_linear_constraints(
    atom_count,
    charge,
    extra_equivalence=None,
    constraint_matrix=None,
    constraint_targets=None,
):
    """Build a deterministic independent ``C q = d`` system."""

    if extra_equivalence is None:
        extra_equivalence = []
    rows = [np.ones(atom_count, dtype=float)]
    targets = [float(charge)]
    labels = ["total_charge"]
    if constraint_matrix is not None:
        matrix = np.asarray(constraint_matrix, dtype=float)
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)
        if constraint_targets is None:
            raise ValueError(
                "RESP constraint targets require a constraint matrix"
            )
        values = np.asarray(constraint_targets, dtype=float).reshape(-1)
        if (
            matrix.ndim != 2
            or matrix.shape[1] != atom_count
            or matrix.shape[0] != values.size
        ):
            raise ValueError(
                "RESP linear constraints have incompatible dimensions"
            )
        if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(values)):
            raise ValueError("RESP linear constraints should be finite")
        rows.extend(matrix)
        targets.extend(values)
        labels.extend(
            f"linear_constraint[{index}]" for index in range(len(values))
        )
    elif constraint_targets is not None:
        raise ValueError("RESP constraint targets require a constraint matrix")

    for group_index, group in enumerate(extra_equivalence):
        indices = [int(index) for index in group]
        if len(indices) < 2 or len(indices) != len(set(indices)):
            raise ValueError(
                f"RESP equivalence group {group_index} should contain "
                "unique atom indices"
            )
        if min(indices) < 0 or max(indices) >= atom_count:
            raise ValueError(
                f"RESP equivalence group {group_index} contains an "
                "out-of-range atom index"
            )
        reference = indices[0]
        for index in indices[1:]:
            row = np.zeros(atom_count, dtype=float)
            row[index] = 1.0
            row[reference] = -1.0
            rows.append(row)
            targets.append(0.0)
            labels.append(f"equivalence[{group_index}]({reference},{index})")

    matrix = np.asarray(rows, dtype=float)
    values = np.asarray(targets, dtype=float)
    scale = max(
        1.0,
        float(np.max(np.abs(matrix))),
        float(np.max(np.abs(values))),
    )
    tolerance = 1.0e-10 * scale
    candidate, _, _, _ = np.linalg.lstsq(matrix, values, rcond=None)
    if float(np.max(np.abs(matrix @ candidate - values))) > tolerance:
        raise ValueError("RESP linear constraints are inconsistent")

    independent_rows = []
    independent_targets = []
    independent_indices = []
    for row_index, (row, target) in enumerate(zip(matrix, values)):
        if independent_rows:
            current = np.asarray(independent_rows, dtype=float)
            coefficients, _, _, _ = np.linalg.lstsq(
                current.T, row, rcond=None
            )
            residual = row - coefficients @ current
            row_scale = max(1.0, float(np.max(np.abs(row))))
            if (
                float(np.max(np.abs(residual)))
                <= 1.0e-10 * row_scale
            ):
                continue
        independent_rows.append(row)
        independent_targets.append(target)
        independent_indices.append(row_index)
    independent = set(independent_indices)
    ledger = [
        {
            "label": label,
            "target": float(target),
            "coefficients": row.tolist(),
            "independent": index in independent,
        }
        for index, (label, row, target) in enumerate(
            zip(labels, matrix, values)
        )
    ]
    rank = len(independent_rows)
    return (
        np.asarray(independent_rows, dtype=float),
        np.asarray(independent_targets, dtype=float),
        {
            "input_constraint_count": int(len(matrix)),
            "constraint_rank": int(rank),
            "dropped_dependent_constraint_count": int(len(matrix) - rank),
            "constraint_ledger": ledger,
        },
    )


def _solve_constrained_quadratic(
    matrix_a, matrix_b, constraints, targets, *, core="python"
):
    matrix_a = np.asarray(matrix_a, dtype=float)
    matrix_b = np.asarray(matrix_b, dtype=float).reshape(-1)
    constraints = np.asarray(constraints, dtype=float)
    targets = np.asarray(targets, dtype=float).reshape(-1)
    atom_count = matrix_b.size
    if matrix_a.shape != (atom_count, atom_count):
        raise ValueError("RESP quadratic matrix has incompatible dimensions")
    if constraints.ndim != 2 or constraints.shape[1] != atom_count:
        raise ValueError("RESP constraint matrix has incompatible dimensions")
    if constraints.shape[0] != targets.size:
        raise ValueError("RESP constraint targets have incompatible dimensions")
    if not all(
        np.all(np.isfinite(value))
        for value in (matrix_a, matrix_b, constraints, targets)
    ):
        raise ValueError("RESP quadratic system should be finite")
    count = constraints.shape[0]
    kkt = np.block(
        [
            [matrix_a, constraints.T],
            [constraints, np.zeros((count, count), dtype=float)],
        ]
    )
    rhs = np.concatenate((matrix_b, targets))
    reference_solution, _, rank, singular_values = np.linalg.lstsq(
        kkt, rhs, rcond=None
    )
    reference_residual = kkt @ reference_solution - rhs
    scale = max(1.0, float(np.max(np.abs(rhs))))
    if (
        rank < kkt.shape[0]
        or float(np.max(np.abs(reference_residual))) > 1.0e-10 * scale
    ):
        raise ValueError(
            "RESP constrained quadratic system is singular or inconsistent"
        )
    if core == "cpp":
        solution = np.asarray(
            _solve_resp_constrained_quadratic_cpp(
                matrix_a.tolist(),
                matrix_b.tolist(),
                constraints.tolist(),
                targets.tolist(),
            ),
            dtype=float,
        )
    elif core == "python":
        solution = reference_solution[:atom_count]
    else:
        raise ValueError("RESP constrained core should be 'python' or 'cpp'")
    charges = solution[:atom_count]
    constraint_residual = constraints @ charges - targets
    condition_number = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and singular_values[-1] > 0
        else None
    )
    return charges, {
        "kkt_rank": int(rank),
        "kkt_size": int(kkt.shape[0]),
        "condition_number": condition_number,
        "max_constraint_residual": float(
            np.max(np.abs(constraint_residual))
        )
        if constraint_residual.size
        else 0.0,
    }


def _fit_constrained_stage(
    assign,
    matrix_a0,
    matrix_b,
    constraints,
    targets,
    initial_charges,
    restraint,
    restrained_atoms,
    core,
):
    charges = np.asarray(initial_charges, dtype=float).reshape(-1)
    restrained_atoms = frozenset(restrained_atoms)
    for iteration in range(1, 1001):
        matrix = np.array(matrix_a0, dtype=float, copy=True)
        for index in restrained_atoms:
            if assign.atoms[index] != "H":
                matrix[index, index] += restraint / np.sqrt(
                    charges[index] * charges[index] + 0.01
                )
        updated, diagnostics = _solve_constrained_quadratic(
            matrix, matrix_b, constraints, targets, core=core
        )
        if float(np.max(np.abs(updated - charges))) <= 1.0e-8:
            diagnostics["iterations"] = iteration
            return updated, diagnostics
        charges = updated
    raise RuntimeError("RESP constrained iteration did not converge")


def _second_stage_groups(assign):
    groups = []
    for atom, element in enumerate(assign.atoms):
        if element != "C":
            continue
        hydrogens = sorted(
            int(neighbor)
            for neighbor in assign.bonds[atom]
            if assign.atoms[int(neighbor)] == "H"
        )
        atom_judge = getattr(assign, "Atom_Judge", None)
        is_c4 = (
            bool(atom_judge(atom, "C4"))
            if callable(atom_judge)
            else len(assign.bonds[atom]) == 4
        )
        is_c3 = (
            bool(atom_judge(atom, "C3"))
            if callable(atom_judge)
            else len(assign.bonds[atom]) == 3
        )
        if is_c4 and hydrogens:
            groups.extend(([atom], hydrogens))
        elif is_c3 and len(hydrogens) == 2:
            groups.extend(([atom], hydrogens))
    return groups


def _assemble_quadratic(
    atom_coordinates_bohr,
    nuclear_charges,
    grid_points_bohr,
    esp_values_au,
):
    coordinates = np.asarray(atom_coordinates_bohr, dtype=float)
    nuclear = np.asarray(nuclear_charges, dtype=float).reshape(-1)
    grids = np.asarray(grid_points_bohr, dtype=float)
    electronic_esp = np.asarray(esp_values_au, dtype=float).reshape(-1)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("RESP atom coordinates should have shape (N, 3)")
    if grids.ndim != 2 or grids.shape[1] != 3:
        raise ValueError("RESP grid points should have shape (M, 3)")
    if len(nuclear) != len(coordinates) or len(electronic_esp) != len(grids):
        raise ValueError("RESP ESP input dimensions do not match")
    distances = np.linalg.norm(
        coordinates[:, None, :] - grids[None, :, :], axis=2
    )
    if np.any(distances <= 0) or not np.all(np.isfinite(distances)):
        raise ValueError("RESP grid points must not coincide with atoms")
    inverse = 1.0 / distances
    matrix_a = inverse @ inverse.T
    molecular_esp = nuclear @ inverse - electronic_esp
    matrix_b = inverse @ molecular_esp
    return matrix_a, matrix_b, molecular_esp, inverse


def _fit_with_constraints(
    assign,
    matrix_a,
    matrix_b,
    charge,
    extra_equivalence,
    constraint_matrix,
    constraint_targets,
    a1,
    a2,
    two_stage,
    only_esp,
    core,
):
    atom_count = len(assign.atoms)
    constraints, targets, diagnostics = _prepare_linear_constraints(
        atom_count,
        charge,
        extra_equivalence=extra_equivalence,
        constraint_matrix=constraint_matrix,
        constraint_targets=constraint_targets,
    )
    charges, initial = _solve_constrained_quadratic(
        matrix_a, matrix_b, constraints, targets, core=core
    )
    diagnostics.update({"initial": initial, "stage1": None, "stage2": None})
    if only_esp:
        diagnostics["max_constraint_residual"] = initial[
            "max_constraint_residual"
        ]
        return charges, diagnostics
    charges, stage1 = _fit_constrained_stage(
        assign,
        matrix_a,
        matrix_b,
        constraints,
        targets,
        charges,
        a1,
        range(atom_count),
        core,
    )
    diagnostics["stage1"] = stage1
    if not two_stage:
        diagnostics["max_constraint_residual"] = stage1[
            "max_constraint_residual"
        ]
        return charges, diagnostics
    groups = _second_stage_groups(assign)
    if not groups:
        diagnostics["max_constraint_residual"] = stage1[
            "max_constraint_residual"
        ]
        return charges, diagnostics
    rows = list(constraints)
    values = list(targets)
    active = set()
    for group in groups:
        active.update(group)
        for index in group[1:]:
            row = np.zeros(atom_count)
            row[index] = 1.0
            row[group[0]] = -1.0
            rows.append(row)
            values.append(0.0)
    for index in range(atom_count):
        if index not in active:
            row = np.zeros(atom_count)
            row[index] = 1.0
            rows.append(row)
            values.append(charges[index])
    stage2_constraints, stage2_targets, stage2_constraint_info = (
        _prepare_linear_constraints(
            atom_count,
            charge,
            constraint_matrix=np.asarray(rows),
            constraint_targets=np.asarray(values),
        )
    )
    charges, stage2 = _fit_constrained_stage(
        assign,
        matrix_a,
        matrix_b,
        stage2_constraints,
        stage2_targets,
        charges,
        a2,
        active,
        core,
    )
    stage2.update(stage2_constraint_info)
    diagnostics["stage2"] = stage2
    diagnostics["max_constraint_residual"] = stage2[
        "max_constraint_residual"
    ]
    return charges, diagnostics


def get_mk_grid(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
    """Compatibility wrapper that routes MK-grid generation to C++."""
    if radius is None:
        radius = {}
    return _generate_resp_mk_grid_cpp(
        list(assign.atoms),
        atom_coordinates_bohr,
        area_density,
        layer,
        radius,
    )


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
    constraint_matrix=None,
    constraint_targets=None,
    return_diagnostics=False,
    core="python",
):
    """Fit RESP charges, using the constrained path when requested."""
    if extra_equivalence is None:
        extra_equivalence = []
    if (
        constraint_matrix is not None
        or constraint_targets is not None
        or return_diagnostics
    ):
        matrix_a, matrix_b, molecular_esp, inverse = _assemble_quadratic(
            atom_coordinates_bohr,
            nuclear_charges,
            grid_points_bohr,
            esp_values_au,
        )
        charges, diagnostics = _fit_with_constraints(
            assign,
            matrix_a,
            matrix_b,
            int(charge),
            extra_equivalence,
            constraint_matrix,
            constraint_targets,
            a1,
            a2,
            two_stage,
            only_esp,
            core,
        )
        if return_diagnostics:
            residual = charges @ inverse - molecular_esp
            rmse = float(np.sqrt(np.mean(residual * residual)))
            reference_rms = float(
                np.sqrt(np.mean(molecular_esp * molecular_esp))
            )
            diagnostics.update(
                {
                    "esp_point_count": int(len(molecular_esp)),
                    "esp_rmse_au": rmse,
                    "esp_relative_rmse": (
                        rmse / reference_rms
                        if reference_rms > np.finfo(float).eps
                        else 0.0
                    ),
                    "esp_mae_au": float(np.mean(np.abs(residual))),
                    "esp_max_abs_error_au": float(
                        np.max(np.abs(residual))
                    ),
                }
            )
            return {
                "charges": charges.tolist(),
                "diagnostics": diagnostics,
            }
        return charges.tolist()
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
    constraint_matrix=None,
    constraint_targets=None,
    core="python",
):
    """Compatibility wrapper that routes RESP debug fitting to C++."""
    if extra_equivalence is None:
        extra_equivalence = []
    if constraint_matrix is not None or constraint_targets is not None:
        result = fit_resp_from_esp(
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
            constraint_matrix=constraint_matrix,
            constraint_targets=constraint_targets,
            return_diagnostics=True,
            core=core,
        )
        return {
            "final_charges": result["charges"],
            "diagnostics": result["diagnostics"],
        }
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
    constraint_matrix=None,
    constraint_targets=None,
    return_diagnostics=False,
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
        constraint_matrix=constraint_matrix,
        constraint_targets=constraint_targets,
        return_diagnostics=return_diagnostics,
        core="cpp",
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
    constraint_matrix=None,
    constraint_targets=None,
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
        constraint_matrix=constraint_matrix,
        constraint_targets=constraint_targets,
        core="cpp",
    )
