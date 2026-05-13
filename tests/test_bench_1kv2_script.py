import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "benchmarks" / "bench_1kv2.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("bench_1kv2", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_bench_1kv2_parser_defaults_include_8a_and_20a():
    script = _load_script()
    args = script.build_parser().parse_args([])

    assert args.repeat == 5
    assert args.n_solvent == 512
    assert args.padding == [8.0, 20.0]
