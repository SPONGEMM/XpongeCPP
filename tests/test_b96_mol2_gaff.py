from pathlib import Path
import subprocess
import textwrap

import pytest
import XpongeCPP as Xponge
from conftest import DATA_1KV2_DIR, original_xponge_repo


DATA_DIR = DATA_1KV2_DIR
XPONGE_REPO = original_xponge_repo()
B96_MOL2 = DATA_DIR / "B96.mol2"
B96_FRCMOD = DATA_DIR / "B96.frcmod"


def _write_xpongecpp_b96(dirname):
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401

    Xponge.load_frcmod(str(B96_FRCMOD))
    mol = Xponge.load_mol2(str(B96_MOL2))
    Xponge.Save_SPONGE_Input(mol, prefix="b96", dirname=str(dirname))
    return mol


def _write_xponge_b96_reference(dirname):
    if not XPONGE_REPO.exists():
        pytest.skip("local Xponge reference repository is not available")
    script = textwrap.dedent(
        f"""
        from pathlib import Path
        import Xponge
        import Xponge.forcefield.amber.gaff  # noqa: F401
        from Xponge.forcefield import amber

        out = Path({str(dirname)!r})
        out.mkdir(parents=True, exist_ok=True)
        amber.load_parameters_from_frcmod({str(B96_FRCMOD)!r}, prefix=False)
        mol = Xponge.load_mol2({str(B96_MOL2)!r})
        Xponge.save_sponge_input(mol, prefix="b96", dirname=str(out))
        """
    )
    result = subprocess.run(
        ["python", "-c", script],
        cwd=XPONGE_REPO,
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        pytest.skip(f"local Xponge reference failed: {result.stderr[-500:]}")


def _header(dirname, key):
    return (dirname / f"b96_{key}.txt").read_text().splitlines()[0]


def _numeric_rows(path, atom_columns, float_columns, unordered_angle=False):
    rows = []
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        words = line.split()
        atoms = tuple(int(words[i]) for i in atom_columns)
        if unordered_angle:
            atoms = (min(atoms[0], atoms[2]), atoms[1], max(atoms[0], atoms[2]))
        floats = tuple(round(float(words[i]), 6) for i in float_columns)
        rows.append((*atoms, *floats))
    return set(rows)


def _dihedral_rows(path):
    rows = []
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        words = line.split()
        atoms = tuple(int(word) for word in words[:4])
        term = (int(words[4]), round(float(words[5]), 6), round(float(words[6]), 6))
        if term in {(2, 10.5, 3.141593), (2, 1.1, 3.141593), (2, 1.0, 3.141593)}:
            atoms = tuple(sorted(atoms))
        rows.append((atoms, term))
    return set(rows)


def test_b96_typed_mol2_and_frcmod_export_expected_headers(tmp_path):
    mol = _write_xpongecpp_b96(tmp_path)

    assert mol.name == "B96"
    assert mol.atom_count == 76
    assert mol.residue_count == 1
    assert mol.residues[0].name == "B"
    assert mol.validate()

    expected_headers = {
        "residue": "76 1",
        "resname": "1",
        "atom_name": "76",
        "atom_type_name": "76",
        "mass": "76",
        "charge": "76",
        "coordinate": "76",
        "LJ": "76 9",
        "bond": "80",
        "angle": "141",
        "dihedral": "221",
        "exclude": "76 403",
        "nb14": "182",
    }
    for key, header in expected_headers.items():
        assert _header(tmp_path, key) == header


def test_b96_sponge_output_matches_xponge_reference(tmp_path):
    ref_dir = tmp_path / "xponge"
    cpp_dir = tmp_path / "xpongecpp"
    _write_xponge_b96_reference(ref_dir)
    _write_xpongecpp_b96(cpp_dir)

    for key in ["residue", "resname", "atom_name", "atom_type_name", "coordinate", "mass", "charge", "LJ"]:
        assert (cpp_dir / f"b96_{key}.txt").read_text().splitlines() == (
            ref_dir / f"b96_{key}.txt"
        ).read_text().splitlines()

    assert _numeric_rows(cpp_dir / "b96_bond.txt", (0, 1), (2, 3)) == _numeric_rows(
        ref_dir / "b96_bond.txt", (0, 1), (2, 3)
    )
    assert _numeric_rows(cpp_dir / "b96_angle.txt", (0, 1, 2), (3, 4), unordered_angle=True) == _numeric_rows(
        ref_dir / "b96_angle.txt", (0, 1, 2), (3, 4), unordered_angle=True
    )
    assert _dihedral_rows(cpp_dir / "b96_dihedral.txt") == _dihedral_rows(ref_dir / "b96_dihedral.txt")
    assert [
        line.split() for line in (cpp_dir / "b96_exclude.txt").read_text().splitlines()
    ] == [line.split() for line in (ref_dir / "b96_exclude.txt").read_text().splitlines()]
    assert _numeric_rows(cpp_dir / "b96_nb14.txt", (0, 1), (2, 3)) == _numeric_rows(
        ref_dir / "b96_nb14.txt", (0, 1), (2, 3)
    )


def test_1kv2_b96_no_solvent_assembly_headers_match_xponge_reference(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401

    Xponge.load_frcmod(str(B96_FRCMOD))
    protein = Xponge.load_pdb(str(DATA_DIR / "1KV2_H.pdb"))
    ligand = Xponge.load_mol2(str(B96_MOL2))

    Xponge.Add_Molecule(protein, ligand)
    Xponge.Save_SPONGE_Input(protein, prefix="b96", dirname=str(tmp_path))

    assert protein.validate()
    assert protein.atom_count == 5869
    assert protein.residue_count == 361
    assert protein.residues[-1].name == "B"
    expected_headers = {
        "residue": "5869 361",
        "resname": "361",
        "atom_name": "5869",
        "atom_type_name": "5869",
        "coordinate": "5869",
        "mass": "5869",
        "charge": "5869",
        "LJ": "5869 15",
        "bond": "5939",
        "angle": "10746",
        "dihedral": "20248",
        "exclude": "5869 32150",
        "nb14": "15465",
    }
    for key, header in expected_headers.items():
        assert _header(tmp_path, key) == header


def test_1kv2_b96_10a_water_ion_headers_match_xponge_reference(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    Xponge.load_frcmod(str(B96_FRCMOD))
    protein = Xponge.load_pdb(str(DATA_DIR / "1KV2_H.pdb"))
    ligand = Xponge.load_mol2(str(B96_MOL2))
    Xponge.Add_Molecule(protein, ligand)
    Xponge.Add_Solvent_Box(
        protein,
        Xponge.get_template_molecule("WAT"),
        10.0,
        tolerance=2.5,
        seed=20260509,
    )
    Xponge.Add_Ions(protein, {"NA": 64, "CL": 52}, seed=20260509)
    Xponge.Save_SPONGE_Input(protein, prefix="b96", dirname=str(tmp_path))

    assert protein.validate()
    assert protein.atom_count == 51963
    assert protein.residue_count == 15803
    assert protein.residue_counts()["B"] == 1
    expected_headers = {
        "residue": "51963 15803",
        "resname": "15803",
        "atom_name": "51963",
        "atom_type_name": "51963",
        "coordinate": "51963",
        "mass": "51963",
        "charge": "51963",
        "LJ": "51963 18",
        "bond": "51917",
        "angle": "10746",
        "dihedral": "20248",
        "exclude": "51963 78128",
        "nb14": "15465",
    }
    for key, header in expected_headers.items():
        assert _header(tmp_path, key) == header


def test_b96_h_mol2_gaff_assign_matches_xponge_atom_types(tmp_path):
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401

    assignment = Xponge.get_assignment_from_mol2(str(DATA_DIR / "B96_H.mol2"), total_charge="sum")
    assignment.determine_atom_type("gaff")

    assert assignment.atom_count == 76
    assert assignment.bond_count == 80
    assert assignment.atom_types == [
        "c", "o", "n", "ca", "ca", "ca", "ca", "ca", "ca", "ca", "ca", "ca", "ca", "n",
        "c2", "c2", "c2", "n2", "na", "ca", "c3", "c3", "c3", "c3", "ca", "ca", "ca",
        "ca", "ca", "c3", "os", "c3", "c3", "n3", "c3", "c3", "os", "c3", "c3", "hn",
        "ha", "ha", "ha", "ha", "ha", "ha", "hn", "ha", "hc", "hc", "hc", "hc", "hc",
        "hc", "hc", "hc", "hc", "ha", "ha", "ha", "ha", "hc", "hc", "hc", "h1", "h1",
        "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1",
    ]

    typed = assignment.to_molecule("B")
    Xponge.Save_Mol2(typed, str(tmp_path / "typed_b96.mol2"))
    reloaded = Xponge.load_mol2(str(tmp_path / "typed_b96.mol2"))
    assert reloaded.atom_count == 76
    assert reloaded.residue_count == 1
    assert [atom.type for atom in reloaded.residues[0].atoms] == assignment.atom_types
