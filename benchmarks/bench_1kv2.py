#!/usr/bin/env python3
"""Smoke and local performance benchmark for the 1KV2 v1 workflow."""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import statistics
import tempfile
import time
from pathlib import Path


def _load_paths_module():
    module_path = Path(__file__).resolve().with_name("_paths.py")
    spec = importlib.util.spec_from_file_location("_xpongecpp_bench_paths", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


PDB_1KV2_H = _load_paths_module().PDB_1KV2_H


def bench_xpongecpp(padding: float, n_solvent: int, n_repeat: int) -> dict[str, float]:
    import XpongeCPP as Xponge

    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    timings: dict[str, list[float]] = {
        "load_pdb": [],
        "solvate": [],
        "add_ions": [],
        "save_sponge_input": [],
        "save_pdb": [],
        "total": [],
    }
    water = Xponge.get_template_molecule("WAT")

    for _ in range(n_repeat):
        start_total = time.perf_counter()

        start = time.perf_counter()
        mol = Xponge.load_pdb(str(PDB_1KV2_H))
        timings["load_pdb"].append(time.perf_counter() - start)

        start = time.perf_counter()
        Xponge.Add_Solvent_Box(mol, water, padding, tolerance=2.5, n_solvent=n_solvent, seed=20260509)
        timings["solvate"].append(time.perf_counter() - start)

        start = time.perf_counter()
        Xponge.Add_Ions(mol, {"NA": min(64, n_solvent // 2), "CL": min(52, n_solvent // 2)})
        timings["add_ions"].append(time.perf_counter() - start)

        with tempfile.TemporaryDirectory() as tmp:
            start = time.perf_counter()
            Xponge.Save_SPONGE_Input(mol, prefix="spg", dirname=tmp)
            timings["save_sponge_input"].append(time.perf_counter() - start)

            start = time.perf_counter()
            Xponge.save_pdb(mol, str(Path(tmp) / "spg.pdb"))
            timings["save_pdb"].append(time.perf_counter() - start)

        timings["total"].append(time.perf_counter() - start_total)

    return {key: statistics.median(values) for key, values in timings.items()}


def try_import_old_xponge():
    try:
        return importlib.import_module("Xponge")
    except Exception:
        return None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--n-solvent", type=int, default=512)
    parser.add_argument("--padding", type=float, nargs="+", default=[8.0, 20.0])
    parser.add_argument("--smoke", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    repeat = 1 if args.smoke else args.repeat
    n_solvent = 64 if args.smoke else args.n_solvent
    paddings = [8.0] if args.smoke else args.padding

    for padding in paddings:
        result = bench_xpongecpp(padding=padding, n_solvent=n_solvent, n_repeat=repeat)
        print(f"XpongeCPP median wall time (seconds), padding={padding:.1f}A")
        for key, value in result.items():
            print(f"{key:18s} {value:.6f}")
        print()

    if try_import_old_xponge() is None:
        print("old_xponge        unavailable; run from an environment with Xponge installed for ratio comparison")


if __name__ == "__main__":
    main()
