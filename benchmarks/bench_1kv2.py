#!/usr/bin/env python3
"""Smoke and local performance benchmark for the 1KV2 v1 workflow."""

from __future__ import annotations

import argparse
import importlib
import statistics
import tempfile
import time
from pathlib import Path


PDB_1KV2_H = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data/1KV2_H.pdb")


def bench_xpongecpp(n_solvent: int, n_repeat: int) -> dict[str, float]:
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
        Xponge.Add_Solvent_Box(mol, water, 8.0, tolerance=2.5, n_solvent=n_solvent)
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--n-solvent", type=int, default=512)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    repeat = 1 if args.smoke else args.repeat
    n_solvent = 64 if args.smoke else args.n_solvent
    result = bench_xpongecpp(n_solvent=n_solvent, n_repeat=repeat)

    print("XpongeCPP median wall time (seconds)")
    for key, value in result.items():
        print(f"{key:18s} {value:.6f}")

    if try_import_old_xponge() is None:
        print("old_xponge        unavailable; run from an environment with Xponge installed for ratio comparison")


if __name__ == "__main__":
    main()
