from pathlib import Path


def test_1kv2_testdata_bundle_exists_in_repo():
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "tests" / "data" / "1kv2"
    assert (data_dir / "1KV2_H.pdb").is_file()
    assert (data_dir / "B96.mol2").is_file()
    assert (data_dir / "B96_H.mol2").is_file()
    assert (data_dir / "B96.frcmod").is_file()
