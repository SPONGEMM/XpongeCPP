import hashlib
import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "benchmarks" / "fetch_external_regression_data.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("fetch_external_regression_data", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rcsb_plan_uses_fixed_official_pdb_download_urls():
    script = _load_script()

    entries = script.planned_rcsb_entries()

    assert [entry.pdb_id for entry in entries] == ["1CRN", "1UBQ", "1AKE", "4HHB", "1BNA"]
    assert entries[0].url == "https://files.rcsb.org/download/1CRN.pdb"
    assert all(entry.cache_name == f"{entry.pdb_id}.pdb" for entry in entries)


def test_cached_file_manifest_records_sha256_and_relative_path(tmp_path):
    script = _load_script()
    payload = b"HEADER    TEST PDB\nEND\n"
    cached = tmp_path / "1CRN.pdb"
    cached.write_bytes(payload)

    entry = script.RcsbEntry("1CRN", "https://files.rcsb.org/download/1CRN.pdb", "1CRN.pdb")
    manifest_entry = script.cached_file_manifest(entry, cached, fetched_at="2026-05-11")

    assert manifest_entry == {
        "pdb_id": "1CRN",
        "url": "https://files.rcsb.org/download/1CRN.pdb",
        "cache_path": "1CRN.pdb",
        "fetched_at": "2026-05-11",
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
    }


def test_chembl_baseline_defaults_target_1000_molecules():
    script = _load_script()

    args = script.build_parser().parse_args([])

    assert args.chembl_limit == 1000
    assert args.chembl_output_dir == REPO_ROOT / "tests" / "data" / "gaff_assign_1000"
    assert args.rcsb_cache_dir == REPO_ROOT / "tests" / "data" / "rcsb_cache"
