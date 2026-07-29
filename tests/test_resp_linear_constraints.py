import numpy as np
import pytest

import XpongeCPP as Xponge
from XpongeCPP.assign import resp as resp_module
from XpongeCPP.assign import resp_core


def _assignment(elements):
    assignment = Xponge.Assign("constraint-test")
    for index, element in enumerate(elements):
        assignment.add_atom(
            element,
            float(index),
            0.0,
            0.0,
            f"{element}{index + 1}",
            0.0,
        )
    return assignment


def _exact_esp_problem(charges):
    charges = np.asarray(charges, dtype=float)
    atom_coordinates = np.asarray(
        [[0.0, 0.0, 0.0], [2.0, 0.2, 0.0], [0.3, 2.1, 0.1]][
            : len(charges)
        ]
    )
    grids = np.asarray(
        [
            [4.0, 0.5, 0.2],
            [-3.0, 1.0, 0.7],
            [1.0, 4.0, -0.5],
            [0.5, -3.0, 1.2],
            [3.5, 3.0, 0.9],
            [-2.0, -2.5, -0.8],
        ]
    )
    nuclear = np.asarray([8.0, 1.0, 1.0][: len(charges)])
    inverse = 1.0 / np.linalg.norm(
        atom_coordinates[:, None, :] - grids[None, :, :], axis=2
    )
    molecular_esp = charges @ inverse
    electronic_esp = nuclear @ inverse - molecular_esp
    return atom_coordinates, nuclear, grids, electronic_esp


def test_constrained_resp_recovers_linear_and_equivalence_targets():
    assignment = _assignment(["O", "H", "H"])
    problem = _exact_esp_problem([-0.4, 0.2, 0.2])

    result = resp_core.fit_resp_from_esp(
        assignment,
        *problem,
        charge=0,
        extra_equivalence=[[1, 2]],
        constraint_matrix=[[1.0, 0.0, 0.0]],
        constraint_targets=[-0.4],
        only_esp=True,
        two_stage=False,
        return_diagnostics=True,
    )

    assert result["charges"] == pytest.approx([-0.4, 0.2, 0.2], abs=1e-12)
    diagnostics = result["diagnostics"]
    assert diagnostics["constraint_rank"] == 3
    assert diagnostics["max_constraint_residual"] <= 1e-12
    assert diagnostics["esp_rmse_au"] <= 1e-12
    assert [row["label"] for row in diagnostics["constraint_ledger"]] == [
        "total_charge",
        "linear_constraint[0]",
        "equivalence[0](1,2)",
    ]


def test_constrained_resp_drops_dependent_rows_deterministically():
    constraints, targets, diagnostics = (
        resp_core._prepare_linear_constraints(
            3,
            0,
            constraint_matrix=[
                [1.0, 1.0, 1.0],
                [1.0, 0.0, 0.0],
            ],
            constraint_targets=[0.0, 0.25],
        )
    )

    assert constraints.shape == (2, 3)
    assert targets == pytest.approx([0.0, 0.25])
    assert diagnostics["input_constraint_count"] == 3
    assert diagnostics["constraint_rank"] == 2
    assert diagnostics["dropped_dependent_constraint_count"] == 1
    assert diagnostics["constraint_ledger"][1]["independent"] is False


def test_constrained_resp_rejects_conflicting_or_nonfinite_constraints():
    with pytest.raises(ValueError, match="constraints are inconsistent"):
        resp_core._prepare_linear_constraints(
            2,
            0,
            constraint_matrix=[[1.0, 0.0], [1.0, 0.0]],
            constraint_targets=[0.0, 1.0],
        )
    with pytest.raises(ValueError, match="should be finite"):
        resp_core._prepare_linear_constraints(
            2,
            0,
            constraint_matrix=[[np.nan, 0.0]],
            constraint_targets=[0.0],
        )


def test_constrained_resp_rejects_singular_quadratic_system():
    with pytest.raises(
        ValueError, match="quadratic system is singular or inconsistent"
    ):
        resp_core._solve_constrained_quadratic(
            np.zeros((2, 2)),
            np.zeros(2),
            np.ones((1, 2)),
            np.zeros(1),
        )


def test_constrained_resp_is_stable_under_atom_permutation():
    target = np.asarray([-0.4, 0.2, 0.2])
    assignment = _assignment(["O", "H", "H"])
    problem = _exact_esp_problem(target)
    baseline = resp_core.fit_resp_from_esp(
        assignment,
        *problem,
        charge=0,
        extra_equivalence=[[1, 2]],
        constraint_matrix=[[1.0, 0.0, 0.0]],
        constraint_targets=[-0.4],
        only_esp=True,
    )

    permutation = np.asarray([2, 0, 1])
    inverse_permutation = np.argsort(permutation)
    permuted_assignment = _assignment(["H", "O", "H"])
    permuted_problem = (
        problem[0][permutation],
        problem[1][permutation],
        problem[2],
        problem[3],
    )
    permuted = resp_core.fit_resp_from_esp(
        permuted_assignment,
        *permuted_problem,
        charge=0,
        extra_equivalence=[[0, 2]],
        constraint_matrix=[[0.0, 1.0, 0.0]],
        constraint_targets=[-0.4],
        only_esp=True,
    )

    assert np.asarray(permuted)[inverse_permutation] == pytest.approx(
        baseline, abs=1e-12
    )


def test_constrained_resp_python_and_cpp_solvers_agree_on_cached_esp():
    assignment = _assignment(["O", "H", "H"])
    problem = _exact_esp_problem([-0.4, 0.2, 0.2])
    kwargs = {
        "charge": 0,
        "extra_equivalence": [[1, 2]],
        "constraint_matrix": [[1.0, 0.0, 0.0]],
        "constraint_targets": [-0.4],
        "only_esp": False,
        "two_stage": False,
        "return_diagnostics": True,
    }

    python_result = resp_core.fit_resp_from_esp(
        assignment, *problem, core="python", **kwargs
    )
    cpp_result = resp_core.fit_resp_from_esp(
        assignment, *problem, core="cpp", **kwargs
    )

    assert cpp_result["charges"] == pytest.approx(
        python_result["charges"], abs=1e-10
    )
    assert cpp_result["diagnostics"]["max_constraint_residual"] <= 1e-10


def test_resp_preflights_constraints_before_qm_backend(monkeypatch):
    assignment = _assignment(["H", "H"])
    called = False

    def forbidden_backend(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("QM backend should not run")

    monkeypatch.setattr(resp_module, "_build_backend_payload", forbidden_backend)
    with pytest.raises(ValueError, match="constraints are inconsistent"):
        resp_module.resp_fit(
            assignment,
            charge=0,
            constraint_matrix=[[1.0, 0.0], [1.0, 0.0]],
            constraint_targets=[0.0, 1.0],
        )
    assert called is False
