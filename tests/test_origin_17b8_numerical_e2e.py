"""Real-SPONGE numerical gates for the 1.7b8 local-patch metal workflows."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import subprocess

import numpy as np
import pytest

import XpongeCPP as Xponge
from XpongeCPP.assign.resp_core import fit_resp_from_esp
from XpongeCPP.io_bundle import (
    convert_bundle_to_legacy,
    convert_legacy_to_bundle,
)
from XpongeCPP.metal_assignment import (
    BaseForceFieldOverlay,
    BondedParameterOverlay,
    apply,
    build_metal_parameter_patch,
    manual_bonded_terms,
)
from XpongeCPP.metal_assignment.artifacts import (
    DerivedModel,
    ModelAtom,
    ModelBond,
    ModelLink,
)
from XpongeCPP.metal_assignment.force_fit import (
    HessianArtifact,
    seminario_bonded_terms,
)

from io_bundle_fixtures import find_sponge_executable
from origin_test_metal_assignment_apply import (
    _ordinary_molecule,
)
from origin_test_metal_assignment_contracts import (
    _bonded_request,
    _overlay_result,
)


def _preassigned_result(request, result):
    return replace(
        result,
        base_overlay=BaseForceFieldOverlay(
            topology_hash=request.topology.topology_hash,
            parameter_source="xponge:preassigned-molecule",
        ),
        provenance={
            **result.provenance,
            "base_assignment": "preassigned_molecule",
        },
        result_hash="",
    ).with_computed_hash()


def _manual_reference_geometry(request):
    payload = {
        "schema_version": 1,
        "graph_revision": request.topology.graph_revision,
        "input_hash": request.topology.input_hash,
        "coordinate_unit": "angstrom",
        "angle_unit": "rad",
        "geometry_source": "frozen_current_geometry",
        "selections": [{
            "selection_id": "site:heme-fe",
            "center_external_id": "heme-fe",
            "bonds": [
                {
                    "edge_id": "fe-na",
                    "center_external_id": "heme-fe",
                    "neighbor_external_id": "heme-na",
                    "equilibrium": 1.9,
                    "unit": "angstrom",
                },
                {
                    "edge_id": "fe-his",
                    "center_external_id": "heme-fe",
                    "neighbor_external_id": "his-ne2",
                    "equilibrium": 2.1,
                    "unit": "angstrom",
                },
            ],
            "angles": [{
                "edge_id1": "fe-na",
                "edge_id2": "fe-his",
                "neighbor1_external_id": "heme-na",
                "center_external_id": "heme-fe",
                "neighbor2_external_id": "his-ne2",
                "equilibrium": np.pi / 2.0,
                "unit": "rad",
            }],
        }],
    }
    import hashlib
    import json

    payload["artifact_hash"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return payload


def _manual_bonded_case():
    request = _bonded_request()
    terms, report = manual_bonded_terms(
        request.topology,
        bond_force_constant=100.0,
        angle_force_constant=20.0,
        reference_geometry_artifact=_manual_reference_geometry(request),
        site_force_constants={
            "site:heme-fe": {
                "bond_force_constant": 111.0,
                "angle_force_constant": 22.0,
            },
        },
    )
    result = _overlay_result(request)
    result = replace(
        result,
        bonded_overlay=BondedParameterOverlay(
            request.topology.topology_hash,
            terms,
            "manual_bonded:explicit_reference_geometry",
        ),
        fit_reports={**result.fit_reports, "bonded": report},
        provenance={
            **result.provenance,
            "force_method": "manual_bonded",
        },
        result_hash="",
    ).with_computed_hash()
    assert all(
        term["source"] == "manual_bonded:explicit_reference_geometry"
        for term in result.bonded_overlay.terms.values()
    )
    return request, _preassigned_result(request, result)


def _seminario_case():
    request = _bonded_request()
    topology = replace(
        request.topology,
        atoms=tuple(
            replace(atom, coordinates=(0.0, 2.1, 0.4))
            if atom.external_id == "his-ne2"
            else atom
            for atom in request.topology.atoms
        ),
        topology_hash="",
    ).with_computed_hash()
    request = replace(
        request,
        topology=topology,
        projection_hash="",
    ).with_computed_hash()
    atoms_by_id = {
        atom.external_id: atom for atom in request.topology.atoms
    }

    def model_atom(model_id, external_id, serial):
        atom = atoms_by_id[external_id]
        return ModelAtom(
            model_id,
            external_id,
            atom.canonical_atom_id,
            serial,
            atom.element,
            atom.coordinates,
            "core",
            atom.canonical_residue_id,
            atom.chemical_component_id,
            None,
            "",
            "",
            "prepared_topology",
            external_id,
        )

    model = DerivedModel(
        external_id="small:heme-fe",
        site_id="site:heme-fe",
        purpose="small",
        coordinate_unit="angstrom",
        atomic_charge_role="absent",
        electronic_state=None,
        atoms=(
            model_atom("m:fe", "heme-fe", 1),
            model_atom("m:na", "heme-na", 2),
            model_atom("m:his", "his-ne2", 3),
        ),
        bonds=(ModelBond(
            "fe-na",
            ("m:fe", "m:na"),
            1.0,
            "coordination",
            "confirmed_input",
        ),),
        links=(ModelLink(
            "fe-his",
            ("m:fe", "m:his"),
            "coordination",
            "struct_conn",
        ),),
        cut_edges=(),
        charge_accounting={},
        mol2_text="",
        model_hash="",
    )
    model = replace(model, model_hash=model.computed_hash())
    hessian = np.eye(len(model.atoms) * 3, dtype=float) * 0.01
    artifact = HessianArtifact(
        model.external_id,
        model.model_hash,
        tuple(atom.model_atom_id for atom in model.atoms),
        tuple(atom.coordinates for atom in model.atoms),
        hessian,
        "deterministic-e2e",
        "1",
    )
    terms, report = seminario_bonded_terms(model, artifact)
    result = _overlay_result(request)
    result = replace(
        result,
        bonded_overlay=BondedParameterOverlay(
            request.topology.topology_hash,
            terms,
            "seminario:deterministic-e2e:1",
        ),
        fit_reports={**result.fit_reports, "bonded": report},
        provenance={**result.provenance, "force_method": "seminario"},
        result_hash="",
    ).with_computed_hash()
    assert {term["kind"] for term in terms.values()} == {"bond", "angle"}
    return request, _preassigned_result(request, result)


def _constrained_resp_case():
    request = _bonded_request()
    topology = request.topology
    assignment = Xponge.Assign("constrained-resp-e2e")
    for atom in topology.atoms:
        assignment.add_atom(
            atom.element,
            *atom.coordinates,
            atom.external_id,
            0.0,
        )
    target = np.asarray([1.2, -0.3, 0.4, 0.7])
    coordinates = np.asarray(
        [atom.coordinates for atom in topology.atoms], dtype=float
    )
    nuclear = np.asarray([26.0, 7.0, 6.0, 7.0])
    grid = np.asarray([
        [5.0, 4.0, 3.0],
        [-4.0, 3.0, 2.0],
        [2.0, -5.0, 4.0],
        [4.0, 2.0, -5.0],
        [-3.0, -4.0, 3.0],
        [6.0, -2.0, 1.0],
        [-5.0, 1.0, -3.0],
        [1.0, 6.0, -2.0],
    ])
    inverse = 1.0 / np.linalg.norm(
        coordinates[:, None, :] - grid[None, :, :], axis=2
    )
    electronic_esp = nuclear @ inverse - target @ inverse
    fitted = fit_resp_from_esp(
        assignment,
        coordinates,
        nuclear,
        grid,
        electronic_esp,
        charge=2,
        constraint_matrix=[
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 1.0, 0.0],
        ],
        constraint_targets=[1.2, 0.1],
        only_esp=True,
        two_stage=False,
        core="cpp",
        return_diagnostics=True,
    )
    charges = np.asarray(fitted["charges"])
    assert charges == pytest.approx(target, abs=1.0e-9)
    assert fitted["diagnostics"]["max_constraint_residual"] <= 1.0e-10

    result = _overlay_result(request)
    result = replace(
        result,
        base_overlay=replace(
            result.base_overlay,
            charges={
                "heme-na": float(charges[1]),
                "heme-c": float(charges[2]),
                "his-ne2": float(charges[3]),
            },
        ),
        metal_overlay=replace(
            result.metal_overlay,
            charges={"heme-fe": float(charges[0])},
        ),
        fit_reports={
            **result.fit_reports,
            "resp": fitted["diagnostics"],
        },
        provenance={**result.provenance, "charge_method": "constrained_resp"},
        result_hash="",
    ).with_computed_hash()
    return request, _preassigned_result(request, result)


def _materialize(case_factory):
    request, result = case_factory()
    patch = build_metal_parameter_patch(
        request,
        result,
        site_ids=("site:heme-fe",),
    )
    molecule, mapping = _ordinary_molecule(patch)
    for link in patch.required_links:
        molecule.add_coordination_bond(
            mapping[link.atom_ids[0]],
            mapping[link.atom_ids[1]],
        )
    molecule = apply(molecule, patch, mapping).molecule
    molecule.set_box_padding(20.0)
    return molecule


def _write_mdin(case_dir: Path, *, steps: int) -> None:
    keys = sorted(
        path.stem.removeprefix("system_")
        for path in case_dir.glob("system_*.txt")
    )
    lines = [
        'mode = "minimization"',
        f"step_limit = {steps}",
        "cutoff = 8.0",
        "print_zeroth_frame = 1",
        "write_mdout_interval = 1",
        "write_information_interval = 1",
        "write_trajectory_interval = 1",
        'frc = "force.dat"',
        *[
            f'{key}_in_file = "system_{key}.txt"'
            for key in keys
        ],
    ]
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def _run_sponge(executable: Path, case_dir: Path, mdin: str):
    completed = subprocess.run(
        [str(executable), "-mdin", mdin],
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert completed.returncode == 0, (
        f"SPONGE failed in {case_dir}\nstdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    headers = lines[0].split()
    frames = [
        {key: float(value) for key, value in zip(headers, row.split())}
        for row in lines[1:]
    ]
    return frames


@pytest.mark.parametrize(
    ("name", "case_factory", "steps"),
    (
        ("constrained-resp", _constrained_resp_case, 0),
        ("manual-bonded", _manual_bonded_case, 0),
        ("seminario", _seminario_case, 2),
    ),
)
def test_origin_17b8_feature_fixtures_preserve_sponge_numerics(
    tmp_path, name, case_factory, steps
):
    executable = find_sponge_executable()
    if executable is None:
        pytest.skip("SPONGE executable is unavailable")

    raw_dir = tmp_path / name / "raw"
    raw_dir.mkdir(parents=True)
    Xponge.save_sponge_input(
        _materialize(case_factory),
        "system",
        raw_dir,
        format="raw",
    )
    _write_mdin(raw_dir, steps=steps)

    converted_root = tmp_path / name / "converted"
    convert_legacy_to_bundle(raw_dir, converted_root)
    roundtrip_dir = tmp_path / name / "roundtrip"
    convert_bundle_to_legacy(converted_root / "bundle", roundtrip_dir)

    raw_frames = _run_sponge(executable, raw_dir, "mdin.spg.toml")
    roundtrip_frames = _run_sponge(
        executable,
        roundtrip_dir,
        "mdin.legacy.spg.toml",
    )
    assert len(roundtrip_frames) == len(raw_frames)
    for expected, actual in zip(raw_frames, roundtrip_frames):
        assert actual.keys() == expected.keys()
        for key in expected:
            assert actual[key] == pytest.approx(
                expected[key], abs=5.0e-4, rel=2.0e-6
            ), (name, key)

    if steps:
        for filename in (
            "mdcrd.dat",
            "force.dat",
            "restart_coordinate.txt",
        ):
            expected = raw_dir / filename
            actual = roundtrip_dir / filename
            assert expected.is_file() and actual.is_file(), filename
            if filename.endswith(".dat"):
                assert np.fromfile(actual, dtype=np.float32) == pytest.approx(
                    np.fromfile(expected, dtype=np.float32),
                    abs=5.0e-4,
                    rel=2.0e-6,
                )
            else:
                expected_values = np.fromstring(
                    expected.read_text(encoding="utf-8"), sep=" "
                )
                actual_values = np.fromstring(
                    actual.read_text(encoding="utf-8"), sep=" "
                )
                assert actual_values == pytest.approx(
                    expected_values, abs=2.0e-4, rel=2.0e-6
                )
