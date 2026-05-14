from pathlib import Path


def test_1kv2_testdata_bundle_exists_in_repo():
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "tests" / "data" / "1kv2"
    assert (data_dir / "1KV2_H.pdb").is_file()
    assert (data_dir / "B96.mol2").is_file()
    assert (data_dir / "B96_H.mol2").is_file()
    assert (data_dir / "B96.frcmod").is_file()


def test_8ryk_testdata_bundle_exists_in_repo():
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "tests" / "data" / "8ryk"
    assert (data_dir / "8RYK_pdbfixer_H_ed.pdb").is_file()
    assert (data_dir / "spg_init.txt").is_file()
    assert (data_dir / "frcmod_init.txt").is_file()
    assert (data_dir / "edit_struct" / "CCS_3.gaff.mol2").is_file()
    assert (data_dir / "frcmod" / "interactive.frcmod").is_file()
    assert (data_dir / "frcmod" / "interactive.payload.json").is_file()
    assert (data_dir / "sponge" / "input_coordinate.txt").is_file()
