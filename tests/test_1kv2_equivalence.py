from pathlib import Path

import pytest
import XpongeCPP as Xponge
from conftest import DATA_1KV2_DIR, optional_1kv2_baseline_dir


PDB_1KV2_H = DATA_1KV2_DIR / "1KV2_H.pdb"


def _load_forcefield():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401


def test_1kv2_pdb_uses_xponge_terminal_residue_names():
    _load_forcefield()

    mol = Xponge.load_pdb(str(PDB_1KV2_H))

    assert mol.atom_count == 5793
    assert mol.residue_count == 360
    assert mol.residues[0].name == "NMET"
    assert mol.residues[-1].name == "CSER"


def test_1kv2_tip3p_10a_auto_water_count_matches_xponge():
    _load_forcefield()
    mol = Xponge.load_pdb(str(PDB_1KV2_H))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, seed=20260509)

    counts = mol.residue_counts()
    assert counts["WAT"] == 15442
    assert mol.residue_count == 15802
    assert mol.atom_count == 52119
    assert mol.validate()


def test_1kv2_10a_water_ion_counts_match_xponge(tmp_path):
    _load_forcefield()
    mol = Xponge.load_pdb(str(PDB_1KV2_H))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, seed=20260509)
    Xponge.Add_Ions(mol, {"NA": 64, "CL": 52}, seed=20260509)

    counts = mol.residue_counts()
    assert counts["WAT"] == 15326
    assert counts["NA"] == 64
    assert counts["CL"] == 52
    assert mol.atom_count == 51887
    assert mol.residue_count == 15802

    outputs = Xponge.Save_SPONGE_Input(mol, prefix="spg", dirname=str(tmp_path))
    assert (tmp_path / "spg_residue.txt").read_text().splitlines()[0] == "51887 15802"
    assert (tmp_path / "spg_coordinate.txt").read_text().splitlines()[0] == "51887"
    assert set(outputs) == {
        "LJ",
        "angle",
        "atom_name",
        "atom_type_name",
        "bond",
        "charge",
        "coordinate",
        "dihedral",
        "exclude",
        "mass",
        "nb14",
        "residue",
        "resname",
    }


def test_1kv2_10a_sponge_file_headers_match_xponge(tmp_path):
    _load_forcefield()
    mol = Xponge.load_pdb(str(PDB_1KV2_H))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, seed=20260509)
    Xponge.Add_Ions(mol, {"NA": 64, "CL": 52}, seed=20260509)
    Xponge.Save_SPONGE_Input(mol, prefix="spg", dirname=str(tmp_path))

    expected_headers = {
        "residue": "51887 15802",
        "resname": "15802",
        "atom_name": "51887",
        "atom_type_name": "51887",
        "mass": "51887",
        "charge": "51887",
        "coordinate": "51887",
        "LJ": "51887 17",
        "bond": "51837",
        "angle": "10605",
        "dihedral": "20027",
        "exclude": "51887 77725",
        "nb14": "15283",
    }
    for key, header in expected_headers.items():
        assert (tmp_path / f"spg_{key}.txt").read_text().splitlines()[0] == header


def test_1kv2_dihedral_xponge_same_force_groups_match_local_baseline(tmp_path):
    baseline_dir = optional_1kv2_baseline_dir()
    baseline = None if baseline_dir is None else baseline_dir / "spg_dihedral.txt"
    if baseline is None or not baseline.exists():
        pytest.skip("local Xponge baseline output is not available")

    _load_forcefield()
    mol = Xponge.load_pdb(str(PDB_1KV2_H))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, seed=20260509)
    Xponge.Add_Ions(mol, {"NA": 64, "CL": 52}, seed=20260509)
    Xponge.Save_SPONGE_Input(mol, prefix="spg", dirname=str(tmp_path))

    def xponge_equivalent_rows(path):
        rows = []
        for line in path.read_text().splitlines()[1:]:
            if line.strip():
                words = line.split()
                atoms = tuple(int(word) for word in words[:4])
                term = (int(words[4]), round(float(words[5]), 5), round(float(words[6]), 5))
                if term in {(2, 10.5, 3.14159), (2, 1.1, 3.14159), (2, 1.0, 3.14159)}:
                    rows.append((tuple(sorted(atoms)), term))
                else:
                    rows.append((atoms, term))
        return set(rows)

    assert xponge_equivalent_rows(tmp_path / "spg_dihedral.txt") == xponge_equivalent_rows(baseline)
