#!/usr/bin/env python3
"""Local comparison and benchmark helpers for RESP backends and RESP cores."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import time
from pathlib import Path

import XpongeCPP as Xponge
from XpongeCPP.assign import resp as resp_module
from XpongeCPP.assign import resp_core


ROOT = Path(__file__).resolve().parents[1]
FORMAMIDE_MOL2 = ROOT / "tests" / "data" / "mokda_resp" / "formamide_resp.mol2"
MOKDA_MOL2 = ROOT / "tests" / "data" / "mokda_resp" / "CRO_1.charge-preview.capped.mol2"


def _load_fixture(name: str):
    if name == "formamide":
        return Xponge.get_assignment_from_mol2(str(FORMAMIDE_MOL2), total_charge="sum"), 0
    if name == "mokda":
        return Xponge.get_assignment_from_mol2(str(MOKDA_MOL2), total_charge="sum"), 0
    raise ValueError(f"unknown fixture: {name}")


def _summary(values):
    return {
        "median": statistics.median(values) if values else 0.0,
        "min": min(values) if values else 0.0,
        "max": max(values) if values else 0.0,
    }


def _charge_stats(charges):
    return {
        "charge_count": len(charges),
        "total_charge": float(sum(charges)),
        "min_charge": float(min(charges)),
        "max_charge": float(max(charges)),
        "charge_vector": [float(x) for x in charges],
    }


def _delta_stats(lhs, rhs):
    diffs = [float(a - b) for a, b in zip(lhs, rhs)]
    abs_diffs = [abs(x) for x in diffs]
    rms = math.sqrt(sum(x * x for x in diffs) / len(diffs)) if diffs else 0.0
    return {
        "max_abs_delta": max(abs_diffs) if abs_diffs else 0.0,
        "rms_delta": rms,
    }


def _backend_module(name: str):
    try:
        return resp_module._BACKEND_MODULES[name]
    except KeyError as exc:
        raise ValueError(f"unknown backend: {name}") from exc


def _core_fit(core: str, assignment, payload, grids, electron_esp, total_charge):
    kwargs = dict(
        atom_coordinates_bohr=payload["atom_coordinates_bohr"],
        nuclear_charges=payload["nuclear_charges"],
        grid_points_bohr=grids,
        esp_values_au=electron_esp,
        charge=total_charge,
        extra_equivalence=[],
        a1=0.0005,
        a2=0.001,
        two_stage=False,
        only_esp=True,
    )
    if core == "python":
        return resp_core.fit_resp_from_esp(assignment, **kwargs)
    if core == "cpp":
        return resp_core.fit_resp_from_esp_cpp(assignment, **kwargs)
    raise ValueError(f"unknown core: {core}")


def _core_fit_debug(core: str, assignment, payload, grids, electron_esp, total_charge):
    kwargs = dict(
        atom_coordinates_bohr=payload["atom_coordinates_bohr"],
        nuclear_charges=payload["nuclear_charges"],
        grid_points_bohr=grids,
        esp_values_au=electron_esp,
        charge=total_charge,
        extra_equivalence=[],
        a1=0.0005,
        a2=0.001,
        two_stage=True,
        only_esp=False,
    )
    if core == "python":
        return resp_core.fit_resp_from_esp_debug(assignment, **kwargs)
    if core == "cpp":
        return resp_core.fit_resp_from_esp_cpp_debug(assignment, **kwargs)
    raise ValueError(f"unknown core: {core}")


def bench_backend_once(fixture: str, backend: str, *, basis="sto-3g", grid_density=1, grid_cell_layer=1):
    backend_module = _backend_module(backend)
    assignment, total_charge = _load_fixture(fixture)
    payload = backend_module.build_backend_payload(
        assignment, basis, total_charge, 0, False, return_timings=True
    )
    grid_points = resp_core.get_mk_grid(
        assignment,
        payload["atom_coordinates_bohr"],
        area_density=grid_density,
        layer=grid_cell_layer,
    )
    start = time.perf_counter()
    electron_esp = backend_module.compute_esp_on_grid(payload, grid_points)
    esp_time = time.perf_counter() - start
    charges = resp_core.fit_resp_from_esp(
        assignment,
        atom_coordinates_bohr=payload["atom_coordinates_bohr"],
        nuclear_charges=payload["nuclear_charges"],
        grid_points_bohr=grid_points,
        esp_values_au=electron_esp,
        charge=total_charge,
        extra_equivalence=[],
        a1=0.0005,
        a2=0.001,
        two_stage=False,
        only_esp=True,
    )
    timings = {
        "build": float(payload["timings"]["build"]),
        "scf": float(payload["timings"]["scf"]),
        "esp": float(esp_time),
        "total": float(payload["timings"]["build"] + payload["timings"]["scf"] + esp_time),
    }
    return timings, charges


def bench_backend(fixture: str, backend: str, repeat: int, warmup: int):
    for _ in range(warmup):
        bench_backend_once(fixture, backend)
    timings = {"build": [], "scf": [], "esp": [], "total": []}
    charges = None
    for _ in range(repeat):
        run_timings, charges = bench_backend_once(fixture, backend)
        for key, value in run_timings.items():
            timings[key].append(value)
    return {
        "backend": backend,
        "fixture": fixture,
        "timings": {key: _summary(values) for key, values in timings.items()},
        "charge_stats": _charge_stats(charges or []),
    }


def compare_backends(fixture: str, lhs: str, rhs: str, repeat: int, warmup: int):
    lhs_result = bench_backend(fixture, lhs, repeat, warmup)
    rhs_result = bench_backend(fixture, rhs, repeat, warmup)
    deltas = _delta_stats(
        lhs_result["charge_stats"]["charge_vector"],
        rhs_result["charge_stats"]["charge_vector"],
    )
    total_ratio = (
        rhs_result["timings"]["total"]["median"] / lhs_result["timings"]["total"]["median"]
        if lhs_result["timings"]["total"]["median"]
        else float("inf")
    )
    return {
        "mode": "backend-compare",
        "fixture": fixture,
        "lhs": lhs_result,
        "rhs": rhs_result,
        "delta_stats": deltas,
        "speed_ratio_rhs_over_lhs": total_ratio,
    }


def bench_core_once(fixture: str, core: str, *, basis="sto-3g", grid_density=1, grid_cell_layer=1):
    assignment, total_charge = _load_fixture(fixture)
    payload = resp_module.pyscf_backend.build_backend_payload(assignment, basis, total_charge, 0, False)
    timings = {"grid": 0.0, "assembly": 0.0, "stage1": 0.0, "stage2": 0.0, "fit": 0.0, "total": 0.0}
    total_start = time.perf_counter()
    start = time.perf_counter()
    if core == "python":
        grids = resp_core.get_mk_grid(
            assignment, payload["atom_coordinates_bohr"], area_density=grid_density, layer=grid_cell_layer
        )
    elif core == "cpp":
        grids = resp_core.get_mk_grid_cpp(
            assignment, payload["atom_coordinates_bohr"], area_density=grid_density, layer=grid_cell_layer
        )
    else:
        raise ValueError(f"unknown core: {core}")
    timings["grid"] = time.perf_counter() - start
    electron_esp = resp_module.pyscf_backend.compute_esp_on_grid(payload, grids)
    debug = _core_fit_debug(core, assignment, payload, grids, electron_esp, total_charge)
    timings["assembly"] = float(debug["timings"]["assembly"])
    timings["stage1"] = float(debug["timings"]["stage1"])
    timings["stage2"] = float(debug["timings"]["stage2"])
    timings["fit"] = timings["assembly"] + timings["stage1"] + timings["stage2"]
    timings["total"] = time.perf_counter() - total_start
    return timings, debug


def bench_core(fixture: str, core: str, repeat: int, warmup: int):
    for _ in range(warmup):
        bench_core_once(fixture, core)
    timings = {"grid": [], "assembly": [], "stage1": [], "stage2": [], "fit": [], "total": []}
    debug = None
    for _ in range(repeat):
        run_timings, debug = bench_core_once(fixture, core)
        for key, value in run_timings.items():
            timings[key].append(value)
    return {
        "core": core,
        "fixture": fixture,
        "timings": {key: _summary(values) for key, values in timings.items()},
        "charge_stats": _charge_stats((debug or {})["final_charges"] if debug else []),
    }


def compare_cores(fixture: str, lhs: str, rhs: str, repeat: int, warmup: int):
    lhs_speed = bench_core(fixture, lhs, repeat, warmup)
    rhs_speed = bench_core(fixture, rhs, repeat, warmup)
    assignment, total_charge = _load_fixture(fixture)
    payload = resp_module.pyscf_backend.build_backend_payload(assignment, "sto-3g", total_charge, 0, False)
    grids_python = resp_core.get_mk_grid(assignment, payload["atom_coordinates_bohr"], area_density=1, layer=1)
    grids_cpp = resp_core.get_mk_grid_cpp(assignment, payload["atom_coordinates_bohr"], area_density=1, layer=1)
    import numpy as np

    if not np.allclose(grids_cpp, grids_python, atol=1e-10):
        raise RuntimeError("Python and C++ MK grids diverged; cannot isolate RESP-core comparison")
    electron_esp = resp_module.pyscf_backend.compute_esp_on_grid(payload, grids_python)
    lhs_debug = _core_fit_debug(lhs, assignment, payload, grids_python, electron_esp, total_charge)
    rhs_debug = _core_fit_debug(rhs, assignment, payload, grids_python, electron_esp, total_charge)
    total_ratio = (
        rhs_speed["timings"]["total"]["median"] / lhs_speed["timings"]["total"]["median"]
        if lhs_speed["timings"]["total"]["median"]
        else float("inf")
    )
    return {
        "mode": "core-compare",
        "fixture": fixture,
        "lhs": lhs_speed,
        "rhs": rhs_speed,
        "grid_point_count": int(len(grids_python)),
        "esp_delta_stats": _delta_stats(lhs_debug["esp_charges"], rhs_debug["esp_charges"]),
        "stage1_delta_stats": _delta_stats(lhs_debug["stage1_charges"], rhs_debug["stage1_charges"]),
        "final_delta_stats": _delta_stats(lhs_debug["final_charges"], rhs_debug["final_charges"]),
        "speed_ratio_rhs_over_lhs": total_ratio,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["backend", "backend-compare", "core", "core-compare"],
        required=True,
    )
    parser.add_argument("--fixture", choices=["formamide", "mokda"], default="formamide")
    parser.add_argument("--backend", choices=["pyscf", "psi4"], default="pyscf")
    parser.add_argument("--backend-rhs", choices=["pyscf", "psi4"], default="psi4")
    parser.add_argument("--core", choices=["python", "cpp"], default="python")
    parser.add_argument("--core-rhs", choices=["python", "cpp"], default="cpp")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    if args.mode == "backend":
        payload = bench_backend(args.fixture, args.backend, args.repeat, args.warmup)
    elif args.mode == "backend-compare":
        payload = compare_backends(args.fixture, args.backend, args.backend_rhs, args.repeat, args.warmup)
    elif args.mode == "core":
        payload = bench_core(args.fixture, args.core, args.repeat, args.warmup)
    else:
        payload = compare_cores(args.fixture, args.core, args.core_rhs, args.repeat, args.warmup)

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(payload)


if __name__ == "__main__":
    main()
