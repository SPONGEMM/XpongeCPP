#!/usr/bin/env python3
"""Benchmark parser, assignment, topology, solvation, and export migration paths."""

from __future__ import annotations

import argparse
import json
import statistics
import tempfile
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = REPO_ROOT / "benchmarks" / "forcefield_migration_bench.json"
DEFAULT_MD = REPO_ROOT / "benchmarks" / "forcefield_migration_bench.md"
PDB_1KV2_H = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data/1KV2_H.pdb")
DEFAULT_MOL2_DIR = REPO_ROOT / "tests" / "data" / "gaff_assign_100" / "inputs"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--mol2-dir", type=Path, default=DEFAULT_MOL2_DIR)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-md", type=Path, default=DEFAULT_MD)
    return parser


def median_time(repeat: int, func) -> float:
    values = []
    for _ in range(repeat):
        start = time.perf_counter()
        func()
        values.append(time.perf_counter() - start)
    return statistics.median(values)


def first_mol2(mol2_dir: Path) -> Path:
    try:
        return next(iter(sorted(mol2_dir.glob("*.mol2"))))
    except StopIteration as exc:
        raise FileNotFoundError(f"no mol2 files found in {mol2_dir}") from exc


def run_benchmark(repeat: int, smoke: bool, mol2_dir: Path) -> dict[str, float]:
    import XpongeCPP as Xponge

    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    sample_mol2 = first_mol2(mol2_dir)
    water = Xponge.get_template_molecule("WAT")
    n_solvent = 64 if smoke else 512
    results: dict[str, float] = {}

    results["load_mol2"] = median_time(repeat, lambda: Xponge.load_mol2(str(sample_mol2)))
    results["gaff_assign_from_mol2"] = median_time(
        repeat,
        lambda: Xponge.get_assignment_from_mol2(str(sample_mol2), total_charge="sum").determine_atom_type("gaff"),
    )

    if PDB_1KV2_H.exists():
        results["load_pdb"] = median_time(repeat, lambda: Xponge.load_pdb(str(PDB_1KV2_H)))

        def solvate():
            mol = Xponge.load_pdb(str(PDB_1KV2_H))
            Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, n_solvent=n_solvent)

        results["solvation"] = median_time(repeat, solvate)

        def export():
            mol = Xponge.load_pdb(str(PDB_1KV2_H))
            Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, n_solvent=n_solvent)
            with tempfile.TemporaryDirectory() as tmpdir:
                Xponge.Save_SPONGE_Input(mol, prefix="bench", dirname=tmpdir)

        results["export_sponge_input"] = median_time(repeat, export)

    return results


def write_outputs(results: dict[str, float], output_json: Path, output_md: Path) -> None:
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    lines = ["| step | median_seconds |", "| --- | ---: |"]
    lines.extend(f"| {key} | {value:.6f} |" for key, value in results.items())
    output_md.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = build_parser().parse_args()
    repeat = 1 if args.smoke else args.repeat
    results = run_benchmark(repeat, args.smoke, args.mol2_dir)
    write_outputs(results, args.output_json, args.output_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
