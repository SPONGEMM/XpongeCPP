#!/usr/bin/env python3
"""Local comparison and benchmark helpers for RESP backends and the C++ RESP core."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import time
from pathlib import Path

import XpongeCPP as Xponge
from XpongeCPP.qm import compute_esp_on_grid as qm_compute_esp_on_grid
from XpongeCPP.qm import run_scf as qm_run_scf
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


def _payload_from_scf_result(scf_result):
    payload = {
        "atom_symbols": scf_result.atom_symbols,
        "atom_coordinates_bohr": scf_result.coordinates_bohr,
        "nuclear_charges": scf_result.nuclear_charges,
        "charge": scf_result.charge,
        "spin": scf_result.spin,
        "scf_result": scf_result,
    }
    if scf_result.timings:
        payload["timings"] = dict(scf_result.timings)
    return payload


def _core_fit_debug(assignment, payload, grids, electron_esp, total_charge):
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
    return resp_core.fit_resp_from_esp_debug(assignment, **kwargs)


def bench_backend_once(fixture: str, backend: str, *, basis="sto-3g", grid_density=1, grid_cell_layer=1):
    assignment, total_charge = _load_fixture(fixture)
    scf_result = qm_run_scf(
        assignment,
        backend=backend,
        basis=basis,
        charge=total_charge,
        spin=0,
        optimize_geometry=False,
        return_timings=True,
    )
    payload = _payload_from_scf_result(scf_result)
    grid_points = resp_core.get_mk_grid(
        assignment,
        payload["atom_coordinates_bohr"],
        area_density=grid_density,
        layer=grid_cell_layer,
    )
    start = time.perf_counter()
    electron_esp = qm_compute_esp_on_grid(scf_result, grid_points).electronic_esp_au
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


def bench_core_once(fixture: str, *, basis="sto-3g", grid_density=1, grid_cell_layer=1):
    assignment, total_charge = _load_fixture(fixture)
    scf_result = qm_run_scf(
        assignment,
        backend="pyscf",
        basis=basis,
        charge=total_charge,
        spin=0,
        optimize_geometry=False,
        return_timings=False,
    )
    payload = _payload_from_scf_result(scf_result)
    timings = {"grid": 0.0, "assembly": 0.0, "stage1": 0.0, "stage2": 0.0, "fit": 0.0, "total": 0.0}
    total_start = time.perf_counter()
    start = time.perf_counter()
    grids = resp_core.get_mk_grid(
        assignment, payload["atom_coordinates_bohr"], area_density=grid_density, layer=grid_cell_layer
    )
    timings["grid"] = time.perf_counter() - start
    electron_esp = qm_compute_esp_on_grid(scf_result, grids).electronic_esp_au
    debug = _core_fit_debug(assignment, payload, grids, electron_esp, total_charge)
    timings["assembly"] = float(debug["timings"]["assembly"])
    timings["stage1"] = float(debug["timings"]["stage1"])
    timings["stage2"] = float(debug["timings"]["stage2"])
    timings["fit"] = timings["assembly"] + timings["stage1"] + timings["stage2"]
    timings["total"] = time.perf_counter() - total_start
    return timings, debug


def bench_core(fixture: str, repeat: int, warmup: int):
    for _ in range(warmup):
        bench_core_once(fixture)
    timings = {"grid": [], "assembly": [], "stage1": [], "stage2": [], "fit": [], "total": []}
    debug = None
    for _ in range(repeat):
        run_timings, debug = bench_core_once(fixture)
        for key, value in run_timings.items():
            timings[key].append(value)
    return {
        "core": "cpp",
        "fixture": fixture,
        "timings": {key: _summary(values) for key, values in timings.items()},
        "charge_stats": _charge_stats((debug or {})["final_charges"] if debug else []),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["backend", "backend-compare", "core"],
        required=True,
    )
    parser.add_argument("--fixture", choices=["formamide", "mokda"], default="formamide")
    parser.add_argument("--backend", choices=["pyscf", "psi4"], default="pyscf")
    parser.add_argument("--backend-rhs", choices=["pyscf", "psi4"], default="psi4")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    if args.mode == "backend":
        payload = bench_backend(args.fixture, args.backend, args.repeat, args.warmup)
    elif args.mode == "backend-compare":
        payload = compare_backends(args.fixture, args.backend, args.backend_rhs, args.repeat, args.warmup)
    else:
        payload = bench_core(args.fixture, args.repeat, args.warmup)

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(payload)


if __name__ == "__main__":
    main()
