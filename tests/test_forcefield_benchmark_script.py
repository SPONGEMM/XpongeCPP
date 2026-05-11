import importlib.util
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "benchmarks" / "bench_forcefield_migration.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("bench_forcefield_migration", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_benchmark_parser_defaults_to_five_repeats_and_json_summary(tmp_path):
    script = _load_script()

    args = script.build_parser().parse_args(["--smoke", "--output-json", str(tmp_path / "bench.json")])

    assert args.repeat == 5
    assert args.smoke is True
    assert args.output_json == tmp_path / "bench.json"


def test_write_benchmark_outputs_json_and_markdown(tmp_path):
    script = _load_script()
    results = {"load_mol2": 0.001, "gaff_assign_from_mol2": 0.002}

    script.write_outputs(results, tmp_path / "bench.json", tmp_path / "bench.md")

    assert json.loads((tmp_path / "bench.json").read_text()) == results
    assert "| load_mol2 | 0.001000 |" in (tmp_path / "bench.md").read_text()
