#!/usr/bin/env python3
"""Provider-neutral 1KV2 assembly benchmark used by release parity gates."""

from __future__ import annotations

import argparse
import json
import statistics
import tempfile
import time
from pathlib import Path

import Xponge
import Xponge.forcefield.amber.ff14sb  # noqa: F401
import Xponge.forcefield.amber.tip3p  # noqa: F401


def _median(values):
    return statistics.median(values)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--n-solvent", type=int, default=0)
    args = parser.parse_args()

    timings = {"load_pdb": [], "solvate": [], "save_sponge_input": [], "total": []}
    atom_counts = []
    water = Xponge.get_template_molecule("WAT")
    for _ in range(args.repeat):
        total_start = time.perf_counter()
        start = time.perf_counter()
        molecule = Xponge.load_pdb(str(args.pdb))
        timings["load_pdb"].append(time.perf_counter() - start)

        start = time.perf_counter()
        solvent_options = {"tolerance": 2.5}
        if args.n_solvent > 0:
            solvent_options["n_solvent"] = args.n_solvent
        Xponge.Add_Solvent_Box(molecule, water, 8.0, **solvent_options)
        timings["solvate"].append(time.perf_counter() - start)

        with tempfile.TemporaryDirectory(prefix="xponge-parity-bench-") as output:
            start = time.perf_counter()
            Xponge.Save_SPONGE_Input(molecule, prefix="input", dirname=output)
            timings["save_sponge_input"].append(time.perf_counter() - start)
        timings["total"].append(time.perf_counter() - total_start)
        atom_counts.append(len(molecule.atoms))

    if len(set(atom_counts)) != 1:
        raise RuntimeError(f"non-deterministic atom counts: {atom_counts}")
    print(json.dumps({
        "implementation": str(getattr(Xponge, "__mokda_backend__", "xponge") or "xponge"),
        "version": str(getattr(Xponge, "__version__", "")),
        "atom_count": atom_counts[0],
        "repeat": args.repeat,
        "median_seconds": {key: _median(values) for key, values in timings.items()},
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
