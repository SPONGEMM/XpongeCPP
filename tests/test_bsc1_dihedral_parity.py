from __future__ import annotations

from io import StringIO

import pytest

import Xponge
import Xponge.forcefield.amber.bsc1  # noqa: F401
import Xponge.forcefield.amber.tip3p  # noqa: F401


MOL2 = """@<TRIPOS>MOLECULE
BSC1_DIHEDRAL
4 3 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 CM 1 TMP 0.0
2 A2 1.0 0.0 0.0 N* 1 TMP 0.0
3 A3 2.0 0.0 0.0 CT 1 TMP 0.0
4 A4 3.0 0.0 0.0 OS 1 TMP 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
@<TRIPOS>SUBSTRUCTURE
1 TMP 1
"""


@pytest.mark.parametrize(
    ("periodicity", "k", "phi0"),
    (
        (1, 1.03, 184.8),
        (2, 1.52, 7.8),
        (3, 0.26, 209.6),
    ),
)
def test_bsc1_cytosine_chi_uses_all_frcmod_terms(
    tmp_path,
    periodicity,
    k,
    phi0,
):
    molecule = Xponge.load_mol2(StringIO(MOL2))
    Xponge.save_sponge_input(molecule, "system", tmp_path, format="raw")
    rows = [
        line.split()
        for line in (tmp_path / "system_dihedral.txt").read_text().splitlines()[1:]
    ]
    assert any(
        int(row[4]) == periodicity
        and float(row[5]) == pytest.approx(k)
        and float(row[6]) == pytest.approx(
            phi0 / 180.0 * Xponge.pi,
            abs=1.0e-6,
        )
        for row in rows
    )


def test_bsc1_completed_cytosine_keeps_chi_frcmod_terms(tmp_path):
    control = Xponge.load_mol2(StringIO(MOL2))
    control_before = tmp_path / "control-before"
    Xponge.save_sponge_input(control, "system", control_before, format="raw")
    control_repeated = tmp_path / "control-repeated"
    Xponge.save_sponge_input(control, "system", control_repeated, format="raw")
    assert (
        control_before / "system_dihedral.txt"
    ).read_text() == (
        control_repeated / "system_dihedral.txt"
    ).read_text()
    pdb = tmp_path / "cytosine.pdb"
    pdb.write_text(
        "\n".join(
            (
                "ATOM      1  O4'  DC A   1       0.000   0.000   0.000  1.00  0.00           O",
                "ATOM      2  C1'  DC A   1       1.000   0.000   0.000  1.00  0.00           C",
                "ATOM      3  N1   DC A   1       2.000   0.000   0.000  1.00  0.00           N",
                "ATOM      4  C6   DC A   1       3.000   0.000   0.000  1.00  0.00           C",
                "TER",
                "END",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    molecule = Xponge.load_pdb(pdb)
    control_loaded = tmp_path / "control-loaded"
    Xponge.save_sponge_input(control, "system", control_loaded, format="raw")
    assert (
        control_before / "system_dihedral.txt"
    ).read_text() == (
        control_loaded / "system_dihedral.txt"
    ).read_text()
    molecule.add_missing_atoms()
    control_after = tmp_path / "control-after"
    Xponge.save_sponge_input(control, "system", control_after, format="raw")
    assert (
        control_before / "system_dihedral.txt"
    ).read_text() == (
        control_after / "system_dihedral.txt"
    ).read_text()
    Xponge.save_sponge_input(molecule, "system", tmp_path, format="raw")
    names = (tmp_path / "system_atom_name.txt").read_text().splitlines()[1:]
    atom_ids = {name: names.index(name) for name in ("O4'", "C1'", "N1", "C6")}
    expected_atoms = [
        atom_ids["O4'"],
        atom_ids["C1'"],
        atom_ids["N1"],
        atom_ids["C6"],
    ]
    matching_rows = [
        line.split()
        for line in (tmp_path / "system_dihedral.txt").read_text().splitlines()[1:]
        if [int(value) for value in line.split()[:4]] == expected_atoms
    ]
    assert [int(row[4]) for row in matching_rows] == [1, 2, 3]
