from pathlib import Path
import subprocess
import textwrap

import pytest
import XpongeCPP as Xponge


DATA_DIR = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data")
XPONGE_REPO = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge")
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
