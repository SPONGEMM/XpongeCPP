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
CHARMM36_ITP = REPO_ROOT / "third_party" / "xponge_reference_forcefield" / "charmm" / "charmm36" / "forcefield.itp"
OPLSAAM_ITP = REPO_ROOT / "third_party" / "xponge_reference_forcefield" / "opls" / "oplsaam" / "forcefield.itp"


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
    results["parse_charmm36_itp"] = median_time(repeat, lambda: Xponge.load_gromacs_topology_file(str(CHARMM36_ITP)))
    results["parse_oplsaam_itp"] = median_time(repeat, lambda: Xponge.load_opls_itp_file(str(OPLSAAM_ITP)))

    pairwise_mol2 = """@<TRIPOS>MOLECULE
PAIRWISE
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A 0.0 0.0 0.0 S 1 SOL 0.0
2 B 1.0 0.0 0.0 S 1 SOL 0.0
@<TRIPOS>BOND
1 1 2 1
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        sw = tmp / "mw.sw"
        sw.write_text(
            "S-S 1.1 2.2 3.3 4.0 5.0 6.6 7.7 8.8 0.0 0.0\n"
            "S-S-S 0.0 0.0 3.3 0.0 0.0 0.0 0.0 0.0 9.9 10.1\n"
        )
        edip = tmp / "si.edip"
        edip.write_text(
            "S-S 1.1 2.2 3.3 4.4 5.5 6.6 7.7 8.8 0.0 0.0 0.0 12.12 0.0 0.0 0.0 0.0 0.0\n"
            "S-S-S 0.0 0.0 0.0 0.0 0.0 0.0 9.9 10.1 11.11 12.12 13.13 0.0 14.14 15.15 16.16 17.17 18.18\n"
        )

        def export_pairwise():
            from io import StringIO

            mol = Xponge.load_mol2(StringIO(pairwise_mol2))
            Xponge.load_sw_parameter_file(str(sw), mol)
            Xponge.load_edip_parameter_file(str(edip), mol)
            Xponge.Save_SPONGE_Input(mol, prefix="pairwise", dirname=tmpdir)

        results["export_sw_edip"] = median_time(repeat, export_pairwise)

        def export_softcore():
            from io import StringIO

            Xponge.register_amber_lj_parameter("A", "A", 0.2, 1.0)
            Xponge.register_amber_lj_parameter("B", "B", 0.3, 1.5)
            mol = Xponge.load_mol2(StringIO(pairwise_mol2.replace(" S ", " A ")))
            for atom in mol.residues[0].atoms:
                atom.lj_type_b = "B"
                atom.subsys = 1
            mol.enable_lj_soft_core()
            Xponge.Save_SPONGE_Input(mol, prefix="softcore", dirname=tmpdir)

        results["export_lj_softcore"] = median_time(repeat, export_softcore)

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
