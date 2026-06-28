import math
import io
import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
import XpongeCPP as Xponge
from conftest import original_xponge_repo
from XpongeCPP.qm import compute_esp_on_grid as qm_compute_esp_on_grid
from XpongeCPP.qm import compute_hessian as qm_compute_hessian
from XpongeCPP.qm import get_backend as qm_get_backend
from XpongeCPP.qm import get_capabilities as qm_get_capabilities
from XpongeCPP.qm import optimize_geometry as qm_optimize_geometry
from XpongeCPP.qm import run_scf as qm_run_scf
from XpongeCPP.assign import resp as resp_module
from XpongeCPP.assign import resp_core
from XpongeCPP.qm.capabilities import QMCapabilitySet
from XpongeCPP.qm.errors import QMCapabilityError
from XpongeCPP.qm.models import ESPResult, OptimizationResult, SCFResult
from XpongeCPP.qm import scheduler as qm_scheduler


XPONGE_REPO = original_xponge_repo()
TEST_DATA_DIR = Path(__file__).resolve().parent / "data"
MOKDA_RESP_PREVIEW_MOL2 = TEST_DATA_DIR / "mokda_resp" / "CRO_1.charge-preview.capped.mol2"
FORMAMIDE_RESP_MOL2 = TEST_DATA_DIR / "mokda_resp" / "formamide_resp.mol2"


def _assignment(name, atoms, bonds, formal_charges=None):
    assignment = Xponge.Assign(name)
    for index, element in enumerate(atoms):
        assignment.add_atom(element, float(index), 0.0, 0.0, f"{element}{index + 1}", 0.0)
    for atom1, atom2, order in bonds:
        assignment.add_bond(atom1, atom2, order)
    if formal_charges:
        for atom, charge in formal_charges.items():
            assignment.set_formal_charge(atom, charge)
    return assignment


def _assert_charges_close(actual, expected, tol=1e-6):
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        assert math.isclose(got, want, abs_tol=tol)


def _run_original_xponge(script):
    if not XPONGE_REPO.exists():
        pytest.skip("local Xponge reference repository is not available")
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=XPONGE_REPO,
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        pytest.skip(f"local Xponge reference failed: {result.stderr[-500:]}")
    return result.stdout


def test_tpacm4_calculate_charge_matches_xponge_methane_and_ethane():
    methane = _assignment(
        "methane",
        ["C", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
    )
    methane.calculate_charge("tpacm4")
    _assert_charges_close(methane.charges, [-0.28116, 0.07029, 0.07029, 0.07029, 0.07029])
    assert math.isclose(sum(methane.charges), 0.0, abs_tol=1e-9)

    ethane = _assignment(
        "ethane",
        ["C", "C", "H", "H", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)],
    )
    ethane.calculate_charge("TPACM4")
    _assert_charges_close(
        ethane.charges,
        [-0.205298, -0.205298, 0.068433, 0.068433, 0.068433, 0.068433, 0.068433, 0.068433],
    )
    assert math.isclose(sum(ethane.charges), 0.0, abs_tol=1e-9)


def test_tpacm4_calculate_charge_uses_aromatic_and_functional_group_context():
    benzene = _assignment(
        "benzene",
        ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 2),
            (4, 5, 1),
            (5, 0, 2),
            (0, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
            (3, 9, 1),
            (4, 10, 1),
            (5, 11, 1),
        ],
    )
    benzene.calculate_charge("tpacm4")
    _assert_charges_close(benzene.charges[:6], [-0.155438] * 6)
    _assert_charges_close(benzene.charges[6:], [0.155438] * 6)

    methyl_acetate = _assignment(
        "methyl_acetate",
        ["C", "C", "O", "O", "C", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (0, 5, 1),
            (0, 6, 1),
            (0, 7, 1),
            (1, 2, 2),
            (1, 3, 1),
            (3, 4, 1),
            (4, 8, 1),
            (4, 9, 1),
            (4, 10, 1),
        ],
    )
    methyl_acetate.calculate_charge("tpacm4")
    _assert_charges_close(
        methyl_acetate.charges,
        [-0.403173, 0.812999, -0.551684, -0.366511, -0.018544, 0.086684, 0.086684, 0.086684, 0.088952, 0.088952, 0.088952],
    )


def test_tpacm4_respects_explicit_total_charge_and_formal_charge_property():
    methyl_ammonium = _assignment(
        "methyl_ammonium",
        ["C", "N", "H", "H", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)],
        {1: 1},
    )
    assert methyl_ammonium.formal_charges == [0, 1, 0, 0, 0, 0, 0, 0]

    methyl_ammonium.calculate_charge("tpacm4")
    _assert_charges_close(
        methyl_ammonium.charges,
        [-0.046932, -0.223509, 0.116603, 0.116603, 0.116603, 0.306877, 0.306877, 0.306877],
    )
    assert math.isclose(sum(methyl_ammonium.charges), 1.0, abs_tol=1e-9)


def test_calculate_charge_reports_optional_dependencies_clearly(monkeypatch):
    import builtins

    assignment = _assignment("methane", ["C", "H", "H", "H", "H"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    original_import = builtins.__import__
    default_qm = qm_scheduler.normalize_backend_name(None)

    def fake_import(name, *args, **kwargs):
        if name.startswith("rdkit") or name.startswith("pyscf") or name.startswith("psi4"):
            raise ImportError(name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="RDKit"):
        assignment.calculate_charge("gasteiger")
    expected = "PySCF" if default_qm == "pyscf" else "Psi4"
    with pytest.raises(ImportError, match=expected):
        assignment.calculate_charge("resp")


def test_resp_defaults_to_platform_backend(monkeypatch):
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    calls = []
    default_backend = qm_scheduler.normalize_backend_name(None)

    class FakeBackend:
        @staticmethod
        def build_backend_payload(assign, basis, charge, spin, opt):
            import numpy as np

            calls.append(("build", basis, charge, spin, opt))
            return {
                "atom_coordinates_bohr": np.zeros((3, 3)),
                "nuclear_charges": np.array([8.0, 1.0, 1.0]),
            }

        @staticmethod
        def compute_esp_on_grid(payload, grids, *, memory_limit=None, chunk_policy="auto", safety_factor=0.8):
            import numpy as np

            calls.append(("esp", len(grids), memory_limit, chunk_policy, safety_factor))
            return np.zeros(len(grids))

    monkeypatch.setitem(resp_module._BACKEND_MODULES, default_backend, FakeBackend)
    monkeypatch.setattr(resp_module.resp_core, "get_mk_grid", lambda *args, **kwargs: __import__("numpy").zeros((2, 3)))
    monkeypatch.setattr(resp_module.resp_core, "fit_resp_from_esp", lambda *args, **kwargs: [0.0, 0.0, 0.0])

    charges = resp_module.resp_fit(assignment, basis="sto-3g", charge=0, grid_density=1, grid_cell_layer=1, only_esp=True)

    assert charges == [0.0, 0.0, 0.0]
    assert calls[0] == ("build", "sto-3g", 0, 0, False)
    assert calls[1] == ("esp", 2, None, "auto", 0.8)


def test_resp_forwards_esp_chunking_options(monkeypatch):
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    calls = []

    class FakeBackend:
        @staticmethod
        def build_backend_payload(assign, basis, charge, spin, opt):
            import numpy as np

            return {
                "atom_coordinates_bohr": np.zeros((3, 3)),
                "nuclear_charges": np.array([8.0, 1.0, 1.0]),
            }

        @staticmethod
        def compute_esp_on_grid(payload, grids, *, memory_limit=None, chunk_policy="auto", safety_factor=0.8):
            import numpy as np

            calls.append((memory_limit, chunk_policy, safety_factor))
            return np.zeros(len(grids))

    monkeypatch.setitem(resp_module._BACKEND_MODULES, "pyscf", FakeBackend)
    monkeypatch.setattr(resp_module.resp_core, "get_mk_grid", lambda *args, **kwargs: __import__("numpy").zeros((2, 3)))
    monkeypatch.setattr(resp_module.resp_core, "fit_resp_from_esp", lambda *args, **kwargs: [0.0, 0.0, 0.0])

    resp_module.resp_fit(
        assignment,
        esp_memory_limit="256MB",
        esp_chunk_policy="grid",
        esp_safety_factor=0.5,
        only_esp=True,
    )

    assert calls == [("256MB", "grid", 0.5)]


def test_resp_rejects_unknown_backend():
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    with pytest.raises(ValueError, match="RESP backend should be one of"):
        resp_module.resp_fit(assignment, backend="unknown")


def test_resp_rejects_unknown_core():
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    with pytest.raises(ValueError, match="RESP core should be one of"):
        resp_module.resp_fit(assignment, core="unknown")


def test_resp_windows_hint_mentions_psi4(monkeypatch):
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])

    class FailingBackend:
        @staticmethod
        def build_backend_payload(assign, basis, charge, spin, opt):
            raise ImportError("PySCF is required for RESP charge calculation")

    monkeypatch.setitem(resp_module._BACKEND_MODULES, "pyscf", FailingBackend)
    monkeypatch.setattr(resp_module.sys, "platform", "win32")

    with pytest.raises(ImportError, match="install Psi4"):
        resp_module.resp_fit(assignment, backend="pyscf")


def test_resp_explicit_psi4_backend_reports_missing_dependency(monkeypatch):
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])

    class FailingBackend:
        @staticmethod
        def build_backend_payload(assign, basis, charge, spin, opt):
            raise ImportError("Psi4 is required for RESP charge calculation")

    monkeypatch.setitem(resp_module._BACKEND_MODULES, "psi4", FailingBackend)

    with pytest.raises(ImportError, match="conda-forge"):
        resp_module.resp_fit(assignment, backend="psi4")


def test_qm_scheduler_windows_psi4_hint_mentions_external_install(monkeypatch):
    monkeypatch.setattr(qm_scheduler.sys, "platform", "win32")

    with pytest.raises(ImportError, match="official Psi4 installer"):
        qm_scheduler.backend_import_or_hint("psi4", ImportError("Psi4 is required"))


def test_qm_scheduler_exposes_known_backends():
    assert qm_get_backend("pyscf").name == "pyscf"
    assert qm_get_backend("psi4").name == "psi4"
    with pytest.raises(ValueError, match="QM backend should be one of"):
        qm_get_backend("unknown")


def test_qm_scheduler_default_backend_matches_platform(monkeypatch):
    monkeypatch.setattr(qm_scheduler.sys, "platform", "linux")
    assert qm_scheduler.normalize_backend_name(None) == "pyscf"
    assert qm_get_backend(None).name == "pyscf"

    monkeypatch.setattr(qm_scheduler.sys, "platform", "win32")
    assert qm_scheduler.normalize_backend_name(None) == "psi4"
    assert qm_get_backend(None).name == "psi4"


def test_qm_scheduler_runs_pyscf_scf_and_esp_smoke():
    assignment = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    scf_result = qm_run_scf(
        assignment,
        backend="pyscf",
        basis="sto-3g",
        charge=0,
        spin=0,
        optimize_geometry=False,
        return_timings=True,
    )
    assert scf_result.backend_name == "pyscf"
    assert scf_result.converged
    assert len(scf_result.atom_symbols) == assignment.atom_count
    assert set(scf_result.timings) == {"build", "scf", "total"}
    esp_result = qm_compute_esp_on_grid(scf_result, scf_result.coordinates_bohr[:2])
    assert len(esp_result.electronic_esp_au) == 2
    assert set(esp_result.timings) == {"esp"}
    assert esp_result.diagnostics["mode"] in {"full", "grid_chunk", "shell_grid_chunk"}


def test_qm_scheduler_normalizes_esp_request_options(monkeypatch):
    import numpy as np

    requests = []

    class FakeESPBackend:
        name = "fakeesp"

        @staticmethod
        def capabilities():
            return QMCapabilitySet(supports_scf=False, supports_esp=True)

        @staticmethod
        def compute_esp(scf_result, request):
            requests.append(request)
            return ESPResult(
                grid_points_bohr=np.asarray(request.grid_points_bohr, dtype=float),
                electronic_esp_au=np.zeros(len(request.grid_points_bohr), dtype=float),
            )

    monkeypatch.setitem(qm_scheduler._BACKENDS, "fakeesp", FakeESPBackend)
    scf_result = SCFResult(
        backend_name="fakeesp",
        total_energy=None,
        converged=True,
        coordinates_bohr=np.zeros((0, 3), dtype=float),
        nuclear_charges=np.zeros(0, dtype=float),
        charge=0,
        spin=0,
        atom_symbols=[],
    )

    result = qm_compute_esp_on_grid(
        scf_result,
        np.zeros((3, 3), dtype=float),
        memory_limit="512MB",
        chunk_policy="grid",
        safety_factor=0.5,
    )

    assert len(result.electronic_esp_au) == 3
    assert requests[0].memory_limit_bytes == 512 * 1024 * 1024
    assert requests[0].chunk_policy == "grid"
    assert requests[0].safety_factor == pytest.approx(0.5)


def test_qm_scheduler_runs_pyscf_esp_grid_chunk_mode_smoke():
    assignment = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    scf_result = qm_run_scf(
        assignment,
        backend="pyscf",
        basis="sto-3g",
        charge=0,
        spin=0,
        optimize_geometry=False,
    )
    per_grid_bytes = scf_result.backend_handle["mol"].nao_nr() ** 2 * 8
    esp_result = qm_compute_esp_on_grid(
        scf_result,
        scf_result.coordinates_bohr[:2],
        memory_limit=per_grid_bytes + 1,
        chunk_policy="grid",
        safety_factor=1.0,
    )

    assert len(esp_result.electronic_esp_au) == 2
    assert esp_result.diagnostics["mode"] == "grid_chunk"
    assert esp_result.diagnostics["grid_chunk_size"] == 1


def test_qm_scheduler_runs_pyscf_esp_dual_chunk_mode_smoke():
    assignment = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    scf_result = qm_run_scf(
        assignment,
        backend="pyscf",
        basis="sto-3g",
        charge=0,
        spin=0,
        optimize_geometry=False,
    )
    per_grid_bytes = scf_result.backend_handle["mol"].nao_nr() ** 2 * 8
    esp_result = qm_compute_esp_on_grid(
        scf_result,
        scf_result.coordinates_bohr[:2],
        memory_limit=max(1, per_grid_bytes - 1),
        chunk_policy="auto",
        safety_factor=1.0,
    )

    assert len(esp_result.electronic_esp_au) == 2
    assert esp_result.diagnostics["mode"] == "shell_grid_chunk"
    assert esp_result.diagnostics["shell_block_count"] >= 1


def test_qm_scheduler_reports_capabilities_for_known_backends():
    pyscf_caps = qm_get_capabilities("pyscf")
    psi4_caps = qm_get_capabilities("psi4")
    assert pyscf_caps.supports_scf and pyscf_caps.supports_esp
    assert pyscf_caps.supports_geometry_optimization
    assert pyscf_caps.supports_hessian
    assert psi4_caps.supports_scf and psi4_caps.supports_esp
    assert psi4_caps.supports_geometry_optimization
    assert not psi4_caps.supports_hessian


def test_qm_scheduler_supports_non_resp_geometry_optimization_flow(monkeypatch):
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    calls = []

    class FakeGeometryBackend:
        name = "fakeopt"

        @staticmethod
        def capabilities():
            return QMCapabilitySet(supports_scf=True, supports_esp=False, supports_geometry_optimization=True)

        @staticmethod
        def optimize_geometry(molecule, options, assign=None, return_timings=False):
            calls.append((options.backend, options.basis, options.optimize_geometry, return_timings))
            if assign is not None:
                assign.set_coordinate(0, 1.0, 2.0, 3.0)
            return OptimizationResult(
                optimized_coordinates_angstrom=[(1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)],
                converged=True,
                final_energy=-1.23,
                timings={"build": 0.01, "scf": 0.02, "total": 0.03},
            )

    monkeypatch.setitem(qm_scheduler._BACKENDS, "fakeopt", FakeGeometryBackend)
    result = qm_optimize_geometry(
        assignment,
        backend="fakeopt",
        basis="sto-3g",
        charge=0,
        spin=0,
        return_timings=True,
    )

    assert result.converged
    assert result.final_energy == pytest.approx(-1.23)
    assert result.optimized_coordinates_angstrom[0] == pytest.approx((1.0, 2.0, 3.0))
    assert result.timings == {"build": 0.01, "scf": 0.02, "total": 0.03}
    assert calls == [("fakeopt", "sto-3g", True, True)]
    assert assignment.coordinates[0] == pytest.approx([1.0, 2.0, 3.0])


def test_qm_scheduler_runs_pyscf_hessian_smoke():
    assignment = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    result = qm_compute_hessian(assignment, backend="pyscf", basis="sto-3g", charge=0, spin=0, return_timings=True)
    assert len(result.atom_symbols) == assignment.atom_count
    assert result.cartesian_hessian_au.shape == (assignment.atom_count, assignment.atom_count, 3, 3)
    assert set(result.timings) == {"build", "scf", "hessian", "total"}


def test_qm_scheduler_rejects_unsupported_hessian_cleanly():
    assignment = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    with pytest.raises(QMCapabilityError, match="does not support Hessian"):
        qm_compute_hessian(assignment, backend="psi4", basis="sto-3g", charge=0, spin=0)


def test_assign_charge_aliases_and_pubchem_signature_are_xponge_compatible(monkeypatch):
    assignment = _assignment("methane", ["C", "H", "H", "H", "H"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assignment.Calculate_Charge("tpacm4")
    _assert_charges_close(assignment.charges, [-0.28116, 0.07029, 0.07029, 0.07029, 0.07029])

    import builtins

    original_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("pubchempy"):
            raise ImportError(name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="PubChemPy"):
        Xponge.get_assignment_from_pubchem("CC", "smiles")


def test_cif_assignment_returns_lattice_info_like_xponge():
    cif = """
data_demo
_cell_length_a    10.0
_cell_length_b    11.0
_cell_length_c    12.0
_cell_angle_alpha 90.0
_cell_angle_beta  91.0
_cell_angle_gamma 92.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_Cartn_x
_atom_site_Cartn_y
_atom_site_Cartn_z
C1 C 0.0 0.0 0.0
H1 H 1.0 0.0 0.0
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif)
    assert assignment.name == "demo"
    assert assignment.atoms == ["C", "H"]
    assert lattice["cell_length"] == [10.0, 11.0, 12.0]
    assert lattice["cell_angle"] == [90.0, 91.0, 92.0]


def test_set_ph_deprotonates_carboxylic_acid_like_xponge_phmodel():
    acetic_acid = _assignment(
        "acetic_acid",
        ["C", "C", "O", "O", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 4, 1), (0, 5, 1), (0, 6, 1), (1, 2, 2), (1, 3, 1), (3, 7, 1)],
    )

    total_charge = acetic_acid.set_ph(7.0)

    assert total_charge == -1
    assert acetic_acid.atoms == ["C", "C", "O", "O", "H", "H", "H"]
    assert acetic_acid.formal_charges == [0, 0, 0, -1, 0, 0, 0]
    assert acetic_acid.bond_count == 6


def test_set_ph_protonates_deprotonated_alcohol_like_xponge_phmodel():
    methoxide = _assignment(
        "methoxide",
        ["C", "O", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        {1: -1},
    )

    total_charge = methoxide.set_ph(7.0)

    assert total_charge == 0
    assert methoxide.atoms == ["C", "O", "H", "H", "H", "H"]
    assert methoxide.formal_charges == [0, 0, 0, 0, 0, 0]
    assert methoxide.bond_count == 5


def test_set_ph_protonates_acetate_like_xponge_phmodel():
    acetate = _assignment(
        "acetate",
        ["C", "C", "O", "O", "H", "H", "H"],
        [(0, 1, 1), (0, 4, 1), (0, 5, 1), (0, 6, 1), (1, 2, 2), (1, 3, 1)],
        {3: -1},
    )

    total_charge = acetate.set_ph(2.0)

    assert total_charge == 0
    assert acetate.atoms == ["C", "C", "O", "O", "H", "H", "H", "H"]
    assert acetate.formal_charges == [0, 0, 0, 0, 0, 0, 0, 0]
    assert acetate.bond_count == 7


def test_set_ph_protonates_phenoxide_like_xponge_phmodel():
    phenoxide = _assignment(
        "phenoxide",
        ["C", "C", "C", "C", "C", "C", "O", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 2),
            (4, 5, 1),
            (5, 0, 2),
            (0, 7, 1),
            (1, 8, 1),
            (2, 9, 1),
            (3, 10, 1),
            (4, 11, 1),
            (5, 6, 1),
        ],
        {6: -1},
    )

    total_charge = phenoxide.set_ph(2.0)

    assert total_charge == 0
    assert phenoxide.atoms == ["C", "C", "C", "C", "C", "C", "O", "H", "H", "H", "H", "H", "H"]
    assert phenoxide.formal_charges == [0] * 13
    assert phenoxide.bond_count == 13


def test_phmodel_typing_matches_original_xponge_reference_for_carboxylic_acid():
    acetic_acid = _assignment(
        "acetic_acid",
        ["C", "C", "O", "O", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 4, 1), (0, 5, 1), (0, 6, 1), (1, 2, 2), (1, 3, 1), (3, 7, 1)],
    )
    current_types = acetic_acid.determine_atom_type("phmodel")
    current = [[index, atom_type] for index, atom_type in enumerate(current_types)]

    script = textwrap.dedent(
        f"""
        import json
        import sys
        import builtins
        sys.path.insert(0, {str(XPONGE_REPO)!r})
        import Xponge

        builtins.set_global_alternative_names = Xponge.set_global_alternative_names
        import Xponge.assign.phmodel

        assign = Xponge.Assign("acetic_acid")
        for index, element in enumerate(["C", "C", "O", "O", "H", "H", "H", "H"]):
            assign.Add_Atom(element, float(index), 0.0, 0.0, f"{{element}}{{index + 1}}", 0.0)
        for atom1, atom2, order in [(0, 1, 1), (0, 4, 1), (0, 5, 1), (0, 6, 1), (1, 2, 2), (1, 3, 1), (3, 7, 1)]:
            assign.add_bond(atom1, atom2, order)
        print(json.dumps(sorted(assign.determine_atom_type("phmodel").items())))
        """
    )
    reference = json.loads(_run_original_xponge(script))

    assert current == reference


def test_determine_connectivity_accepts_xponge_tolerance_signature():
    assignment = Xponge.Assign("water")
    assignment.add_atom("O", 0.0, 0.0, 0.0)
    assignment.add_atom("H", 0.96, 0.0, 0.0)
    assignment.add_atom("H", -0.24, 0.93, 0.0)

    assignment.determine_connectivity(tolerance=1.0)

    assert assignment.bond_count == 2
    assert sorted((min(i, j), max(i, j), order) for i, bonds in enumerate(assignment.bonds) for j, order in bonds.items() if i < j) == [
        (0, 1, -1),
        (0, 2, -1),
    ]


def test_pdb_and_xyz_assignment_entrypoints_accept_xponge_parameters():
    pdb = """\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000                       C
HETATM    2  H1  LIG A   1       1.090   0.000   0.000                       H
HETATM    3  C2  SOL A   2       5.000   0.000   0.000                       C
HETATM    4  H2  SOL A   2       6.090   0.000   0.000                       H
END
"""
    assignment = Xponge.get_assignment_from_pdb(io.StringIO(pdb), only_residue="LIG", bond_tolerance=1.0, total_charge=0)
    assert assignment.name == "LIG"
    assert assignment.atoms == ["C", "H"]
    assert assignment.bond_count == 1

    xyz = """2
fragment
C 0.000 0.000 0.000
H 1.090 0.000 0.000
"""
    assignment = Xponge.get_assignment_from_xyz(io.StringIO(xyz), bond_tolerance=1.0, total_charge=0)
    assert assignment.name == "fragment"
    assert assignment.bond_count == 1


def test_cif_assignment_supports_fractional_coordinates_and_symmetry_basis():
    cif = """
data_frac
_cell_length_a    10.0
_cell_length_b    20.0
_cell_length_c    30.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
_symmetry_equiv_pos_as_xyz
'x,y,z
x+1/2,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.5 0.25 0.1
H1 H 0.6 0.25 0.1
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif, keep_cell_angle=False)
    assert assignment.name == "frac"
    assert assignment.coordinates[0] == pytest.approx([5.0, 5.0, 3.0])


def test_cif_assignment_handles_nonorthogonal_cell_and_multiple_symmetry_ops():
    cif = """
data_skew
_cell_length_a    10.0
_cell_length_b    11.0
_cell_length_c    12.0
_cell_angle_alpha 80.0
_cell_angle_beta  90.0
_cell_angle_gamma 100.0
_symmetry_equiv_pos_as_xyz
'x,y,z
x+1/2,y+1/2,z
x,y,z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.1 0.2 0.3
H1 H 0.2 0.2 0.3
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif)

    assert assignment.name == "skew"
    assert assignment.atoms == ["C", "H"]
    assert assignment.bond_count == 1
    assert lattice["cell_length"] == [10.0, 11.0, 12.0]
    assert lattice["cell_angle"] == [80.0, 90.0, 100.0]
    assert lattice["basis_position"] == [[1, 1, 1], [1.5, 1.5, 1], [1, 1, 1.5]]
    assert assignment.coordinates[0][0] != pytest.approx(1.0)


def test_cif_assignment_handles_richer_fractional_crystal_case():
    cif = """
data_rich
_cell_length_a    9.1
_cell_length_b    10.2
_cell_length_c    11.3
_cell_angle_alpha 75.0
_cell_angle_beta  88.0
_cell_angle_gamma 96.0
_symmetry_equiv_pos_as_xyz
'x,y,z
-x,y+1/2,z+1/2
x+1/2,-y,z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.1 0.2 0.3
O1 O 0.2 0.2 0.3
H1 H 0.25 0.28 0.3
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 O1 1.2
O1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif)

    assert assignment.name == "rich"
    assert assignment.atoms == ["C", "O", "H"]
    assert assignment.bond_count == 2
    assert lattice["cell_length"] == [9.1, 10.2, 11.3]
    assert lattice["cell_angle"] == [75.0, 88.0, 96.0]
    assert lattice["basis_position"] == [[1, 1, 1], [-1, 1.5, 1.5], [1.5, -1, 1.5]]


def test_cif_fractional_assignment_matches_original_xponge_reference():
    cif = """
data_frac
_cell_length_a    10.0
_cell_length_b    20.0
_cell_length_c    30.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
_symmetry_equiv_pos_as_xyz
'x,y,z
x+1/2,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.5 0.25 0.1
H1 H 0.6 0.25 0.1
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif, keep_cell_angle=False)
    current = {
        "name": assignment.name,
        "atoms": assignment.atoms,
        "coordinates": assignment.coordinates,
        "cell_length": lattice["cell_length"],
        "cell_angle": lattice["cell_angle"],
    }

    script = textwrap.dedent(
        f"""
        import json
        import sys
        import tempfile
        from pathlib import Path

        sys.path.insert(0, {str(XPONGE_REPO)!r})
        import Xponge

        cif = {cif!r}
        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = Path(tmpdir) / "frac.cif"
            cif_path.write_text(cif)
            assignment, lattice = Xponge.get_assignment_from_cif(str(cif_path), keep_cell_angle=False)
        print(json.dumps({{
            "name": assignment.name,
            "atoms": assignment.atoms,
            "coordinates": assignment.coordinate.tolist(),
            "cell_length": lattice["cell_length"],
            "cell_angle": lattice["cell_angle"],
        }}, sort_keys=True))
        """
    )
    reference = json.loads(_run_original_xponge(script))

    assert current == reference
    assert assignment.bond_count == 1
    assert lattice["cell_angle"] == [90, 90, 90]
    assert lattice["basis_position"] == [[1, 1, 1], [1.5, 1, 1]]


def test_phmodel_atom_type_and_set_ph_cover_phenol_and_alcohol():
    phenol = _assignment(
        "phenol",
        ["C", "C", "C", "C", "C", "C", "O", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 2),
            (4, 5, 1),
            (5, 0, 2),
            (0, 7, 1),
            (1, 8, 1),
            (2, 9, 1),
            (3, 10, 1),
            (4, 11, 1),
            (5, 6, 1),
            (6, 12, 1),
        ],
    )
    types = phenol.determine_atom_type("phmodel")
    assert types[12] == "A-phenol"
    assert phenol.set_ph(11.0) == -1
    assert phenol.atoms.count("H") == 5

    alcohol = _assignment("methanol", ["C", "O", "H", "H", "H", "H"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1)])
    assert alcohol.determine_atom_type("phmodel")[5] == "A-alcohol"


def test_save_as_mol2_atomtype_argument_and_equal_atoms_api(tmp_path):
    ethane = _assignment(
        "ethane",
        ["C", "C", "H", "H", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)],
    )
    ethane.calculate_charge("tpacm4")
    mol2 = tmp_path / "ethane.mol2"
    ethane.save_as_mol2(str(mol2), atomtype="sybyl")
    text = mol2.read_text()
    assert "@<TRIPOS>UNITY_ATOM_ATTR" not in text

    with pytest.raises(ValueError, match="atomtype"):
        ethane.save_as_mol2(str(tmp_path / "bad.mol2"), atomtype="gaff")

    try:
        groups = ethane.determine_equal_atoms()
    except ImportError as exc:
        pytest.skip(str(exc))
    assert any(set(group) == {0, 1} for group in groups)


def test_resp_uses_pyscf_backend_or_reports_missing_dependency():
    water = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    try:
        water.calculate_charge("resp", basis="sto-3g", charge=0, grid_density=1, grid_cell_layer=1, only_esp=True)
    except ImportError as exc:
        assert "PySCF" in str(exc)
        return
    assert len(water.charges) == 3
    assert math.isclose(sum(water.charges), 0.0, abs_tol=1e-5)


def test_resp_supports_nontrivial_small_molecule_under_pyscf():
    methyl_acetate = _assignment(
        "methyl_acetate",
        ["C", "C", "O", "O", "C", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (0, 5, 1),
            (0, 6, 1),
            (0, 7, 1),
            (1, 2, 2),
            (1, 3, 1),
            (3, 4, 1),
            (4, 8, 1),
            (4, 9, 1),
            (4, 10, 1),
        ],
    )
    try:
        methyl_acetate.calculate_charge(
            "resp",
            basis="sto-3g",
            charge=0,
            grid_density=1,
            grid_cell_layer=1,
            only_esp=True,
            two_stage=False,
        )
    except ImportError as exc:
        assert "PySCF" in str(exc)
        return

    assert len(methyl_acetate.charges) == 11
    assert math.isclose(sum(methyl_acetate.charges), 0.0, abs_tol=1e-5)
    assert methyl_acetate.charges[2] < 0.0
    assert methyl_acetate.charges[3] > 0.0
    assert methyl_acetate.charges[4] < 0.0


def test_resp_supports_formamide_under_pyscf():
    formamide = _assignment(
        "formamide",
        ["C", "O", "N", "H", "H", "H"],
        [(0, 1, 2), (0, 2, 1), (0, 3, 1), (2, 4, 1), (2, 5, 1)],
    )
    try:
        formamide.calculate_charge(
            "resp",
            basis="sto-3g",
            charge=0,
            grid_density=1,
            grid_cell_layer=1,
            only_esp=True,
            two_stage=False,
        )
    except ImportError as exc:
        assert "PySCF" in str(exc)
        return

    assert len(formamide.charges) == 6
    assert math.isclose(sum(formamide.charges), 0.0, abs_tol=1e-5)
    assert formamide.charges[0] < 0.0
    assert formamide.charges[1] > 0.0
    assert formamide.charges[2] < 0.0


def test_resp_supports_formamide_under_psi4_if_available():
    pytest.importorskip("psi4")
    formamide = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    formamide.calculate_charge(
        "resp",
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="psi4",
    )

    assert len(formamide.charges) == formamide.atom_count == 6
    assert math.isclose(sum(formamide.charges), 0.0, abs_tol=1e-5)
    assert min(formamide.charges) < -0.3
    assert max(formamide.charges) > 0.2


def test_resp_psi4_matches_pyscf_on_formamide_if_available():
    pytest.importorskip("psi4")
    formamide_pyscf = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    formamide_psi4 = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    formamide_pyscf.calculate_charge(
        "resp",
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="pyscf",
    )
    formamide_psi4.calculate_charge(
        "resp",
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="psi4",
    )

    assert len(formamide_pyscf.charges) == len(formamide_psi4.charges) == 6
    diffs = [abs(a - b) for a, b in zip(formamide_pyscf.charges, formamide_psi4.charges)]
    assert math.isclose(sum(formamide_psi4.charges), 0.0, abs_tol=1e-5)
    assert max(diffs) <= 5e-3


def test_resp_supports_mokda_preview_mol2_under_psi4_if_available():
    pytest.importorskip("psi4")
    assignment = Xponge.get_assignment_from_mol2(str(MOKDA_RESP_PREVIEW_MOL2), total_charge="sum")
    assignment.calculate_charge(
        "resp",
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="psi4",
    )
    assert len(assignment.charges) == assignment.atom_count
    assert math.isclose(sum(assignment.charges), 0.0, abs_tol=1e-4)
    assert min(assignment.charges) < -0.3
    assert max(assignment.charges) > 0.2


def test_resp_psi4_matches_pyscf_on_mokda_preview_if_available():
    pytest.importorskip("psi4")
    assignment_pyscf = Xponge.get_assignment_from_mol2(str(MOKDA_RESP_PREVIEW_MOL2), total_charge="sum")
    assignment_psi4 = Xponge.get_assignment_from_mol2(str(MOKDA_RESP_PREVIEW_MOL2), total_charge="sum")
    assignment_pyscf.calculate_charge(
        "resp",
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="pyscf",
    )
    assignment_psi4.calculate_charge(
        "resp",
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="psi4",
    )
    diffs = [abs(a - b) for a, b in zip(assignment_pyscf.charges, assignment_psi4.charges)]
    assert math.isclose(sum(assignment_psi4.charges), 0.0, abs_tol=1e-4)
    assert max(diffs) <= 5e-3


def test_resp_orchestration_accepts_legacy_python_core_alias_and_cpp_core():
    formamide_python_alias = _assignment(
        "formamide",
        ["C", "O", "N", "H", "H", "H"],
        [(0, 1, 2), (0, 2, 1), (0, 3, 1), (2, 4, 1), (2, 5, 1)],
    )
    formamide_cpp = _assignment(
        "formamide",
        ["C", "O", "N", "H", "H", "H"],
        [(0, 1, 2), (0, 2, 1), (0, 3, 1), (2, 4, 1), (2, 5, 1)],
    )
    charges_python_alias = resp_module.resp_fit(
        formamide_python_alias,
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="pyscf",
        core="python",
    )
    charges_cpp = resp_module.resp_fit(
        formamide_cpp,
        basis="sto-3g",
        charge=0,
        grid_density=1,
        grid_cell_layer=1,
        only_esp=True,
        two_stage=False,
        backend="pyscf",
        core="cpp",
    )
    _assert_charges_close(charges_cpp, charges_python_alias, tol=1e-8)


def test_resp_debug_view_matches_plain_core():
    formamide = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    payload = resp_module.pyscf_backend.build_backend_payload(formamide, "sto-3g", 0, 0, False)
    grids = resp_core.get_mk_grid(formamide, payload["atom_coordinates_bohr"], area_density=1, layer=1)
    electron_esp = resp_module.pyscf_backend.compute_esp_on_grid(payload, grids)
    plain = resp_core.fit_resp_from_esp(
        formamide,
        atom_coordinates_bohr=payload["atom_coordinates_bohr"],
        nuclear_charges=payload["nuclear_charges"],
        grid_points_bohr=grids,
        esp_values_au=electron_esp,
        charge=0,
        extra_equivalence=[],
        a1=0.0005,
        a2=0.001,
        two_stage=True,
        only_esp=False,
    )
    debug = resp_core.fit_resp_from_esp_debug(
        formamide,
        atom_coordinates_bohr=payload["atom_coordinates_bohr"],
        nuclear_charges=payload["nuclear_charges"],
        grid_points_bohr=grids,
        esp_values_au=electron_esp,
        charge=0,
        extra_equivalence=[],
        a1=0.0005,
        a2=0.001,
        two_stage=True,
        only_esp=False,
    )
    _assert_charges_close(debug["final_charges"], plain, tol=1e-8)
    assert set(debug["timings"]) == {"assembly", "stage1", "stage2", "total"}


def test_resp_second_stage_restrains_grouped_heavy_atoms_only():
    cyclohexane = _assignment(
        "cyclohexane",
        ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (1, 2, 1),
            (2, 3, 1),
            (3, 4, 1),
            (4, 5, 1),
            (5, 0, 1),
            (0, 6, 1),
            (0, 7, 1),
            (1, 8, 1),
            (1, 9, 1),
            (2, 10, 1),
            (2, 11, 1),
            (3, 12, 1),
            (3, 13, 1),
            (4, 14, 1),
            (4, 15, 1),
            (5, 16, 1),
            (5, 17, 1),
        ],
    )
    coordinates_angstrom = [
        (1.214, 0.701, 0.0),
        (0.0, 1.402, 0.0),
        (-1.214, 0.701, 0.0),
        (-1.214, -0.701, 0.0),
        (0.0, -1.402, 0.0),
        (1.214, -0.701, 0.0),
        (2.157, 1.245, 0.0),
        (1.214, 0.701, 1.09),
        (0.0, 2.49, 0.0),
        (0.0, 1.402, 1.09),
        (-2.157, 1.245, 0.0),
        (-1.214, 0.701, 1.09),
        (-2.157, -1.245, 0.0),
        (-1.214, -0.701, 1.09),
        (0.0, -2.49, 0.0),
        (0.0, -1.402, 1.09),
        (2.157, -1.245, 0.0),
        (1.214, -0.701, 1.09),
    ]
    bohr_per_angstrom = 1.0 / 0.529177210903
    atom_coordinates_bohr = [
        tuple(component * bohr_per_angstrom for component in coordinate) for coordinate in coordinates_angstrom
    ]
    grids = resp_core.get_mk_grid(cyclohexane, atom_coordinates_bohr, area_density=0.25, layer=1)
    debug = resp_core.fit_resp_from_esp_debug(
        cyclohexane,
        atom_coordinates_bohr=atom_coordinates_bohr,
        nuclear_charges=[6.0] * 6 + [1.0] * 12,
        grid_points_bohr=grids,
        esp_values_au=[0.0] * len(grids),
        charge=0,
        a1=0.0005,
        a2=0.001,
        two_stage=True,
        only_esp=False,
    )

    assert debug["stage2_restrained_groups"] == [0, 2, 4, 6, 8, 10]
    assert math.isclose(sum(debug["final_charges"]), 0.0, abs_tol=1e-8)


def test_resp_supports_real_mol2_fixture_and_matches_original_xponge():
    assignment = Xponge.get_assignment_from_mol2(str(FORMAMIDE_RESP_MOL2), total_charge="sum")
    assert assignment.atom_count == 6

    try:
        assignment.calculate_charge(
            "resp",
            basis="sto-3g",
            charge=0,
            grid_density=1,
            grid_cell_layer=1,
            only_esp=True,
            two_stage=False,
        )
    except ImportError as exc:
        assert "PySCF" in str(exc)
        return

    assert len(assignment.charges) == assignment.atom_count
    assert math.isclose(sum(assignment.charges), 0.0, abs_tol=1e-4)
    assert min(assignment.charges) < -0.3
    assert max(assignment.charges) > 0.2

    script = f"""
import json
import sys
sys.path.insert(0, {str(XPONGE_REPO)!r})
import Xponge

assign = Xponge.get_assignment_from_mol2({str(FORMAMIDE_RESP_MOL2)!r}, total_charge="sum")
assign.calculate_charge(
    "resp",
    basis="sto-3g",
    charge=0,
    grid_density=1,
    grid_cell_layer=1,
    only_esp=True,
    two_stage=False,
)
print(json.dumps(assign.charge.tolist()))
"""
    stdout = _run_original_xponge(script)
    expected = json.loads(stdout)
    _assert_charges_close(assignment.charges, expected, tol=5e-5)


def test_resp_supports_mokda_preview_mol2_and_matches_original_xponge():
    assignment = Xponge.get_assignment_from_mol2(str(MOKDA_RESP_PREVIEW_MOL2), total_charge="sum")
    assert assignment.atom_count == 49

    try:
        assignment.calculate_charge(
            "resp",
            basis="sto-3g",
            charge=0,
            grid_density=1,
            grid_cell_layer=1,
            only_esp=True,
            two_stage=False,
        )
    except ImportError as exc:
        assert "PySCF" in str(exc)
        return

    assert len(assignment.charges) == assignment.atom_count
    assert math.isclose(sum(assignment.charges), 0.0, abs_tol=1e-4)
    assert min(assignment.charges) < -0.3
    assert max(assignment.charges) > 0.2

    script = f"""
import json
import sys
sys.path.insert(0, {str(XPONGE_REPO)!r})
import Xponge

assign = Xponge.get_assignment_from_mol2({str(MOKDA_RESP_PREVIEW_MOL2)!r}, total_charge="sum")
assign.calculate_charge(
    "resp",
    basis="sto-3g",
    charge=0,
    grid_density=1,
    grid_cell_layer=1,
    only_esp=True,
    two_stage=False,
)
print(json.dumps(assign.charge.tolist()))
"""
    stdout = _run_original_xponge(script)
    expected = json.loads(stdout)
    _assert_charges_close(assignment.charges, expected, tol=5e-5)



def test_assign_rule_custom_registry_supports_xponge_pure_string_semantics():
    rule = Xponge.AssignRule("unit_custom_rule", pure_string=True)
    calls = []

    def pre_action(assign):
        calls.append(("pre", assign.name))

    def post_action(assign):
        calls.append(("post", assign.name))

    rule.set_pre_action(pre_action)
    rule.set_post_action(post_action)

    @rule.add_rule("hydrogen", priority=10)
    def is_hydrogen(atom, assign):
        return assign.atoms[atom] == "H"

    @rule.add_rule("heavy", priority=0)
    def is_heavy(atom, assign):
        return True

    methane = _assignment(
        "methane",
        ["C", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
    )

    atom_types = methane.determine_atom_type("unit_custom_rule")

    assert atom_types == ["heavy", "hydrogen", "hydrogen", "hydrogen", "hydrogen"]
    assert methane.atom_types == ["", "", "", "", ""]
    assert calls == [("pre", "methane"), ("post", "methane")]


def test_assign_rule_custom_registry_can_assign_in_place():
    rule = Xponge.AssignRule("unit_inplace_rule")

    @rule.add_rule("h_custom", priority=10)
    def is_hydrogen(atom, assign):
        return assign.atoms[atom] == "H"

    @rule.add_rule("x_custom")
    def fallback(atom, assign):
        return True

    water = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])

    assert water.determine_atom_type(rule) is None
    assert water.atom_types == ["x_custom", "h_custom", "h_custom"]


def test_determine_bond_order_accepts_custom_penalty_scores():
    hydroxide = _assignment("hydroxide", ["O", "H"], [(0, 1, -1)])

    success = hydroxide.determine_bond_order(
        penalty_scores=[{1: 0}, {1: 0}],
        total_charge=-1,
    )

    assert success
    assert hydroxide.bonds[0][1] == 1
    assert hydroxide.formal_charges == [-1, 0]


def test_determine_bond_order_extra_criteria_filters_candidates():
    water = _assignment("water", ["O", "H", "H"], [(0, 1, -1), (0, 2, -1)])
    seen = []

    def reject_all(assign):
        seen.append(sorted((min(i, j), max(i, j), order) for i, bonds in enumerate(assign.bonds) for j, order in bonds.items() if i < j))
        return False

    assert not water.determine_bond_order(extra_criteria=reject_all)
    assert seen

    water = _assignment("water", ["O", "H", "H"], [(0, 1, -1), (0, 2, -1)])
    assert water.determine_bond_order(extra_criteria=lambda assign: all(order == 1 for order in assign.bonds[0].values()))
    assert water.bonds[0][1] == 1
    assert water.bonds[0][2] == 1
