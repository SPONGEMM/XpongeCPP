from io import StringIO

import pytest

import XpongeCPP as Xponge


def test_empty_molecule_rejects_box_padding():
    mol = Xponge.Molecule("EMPTY")
    with pytest.raises(ValueError, match="at least one atom"):
        Xponge.Set_Box_Padding(mol, 1.0)


def test_bad_pdb_coordinate_line_reports_parse_error():
    bad = "ATOM      1  N   ALA A   1       BAD     0.000   0.000  1.00  0.00           N\n"
    with pytest.raises(ValueError):
        Xponge.load_pdb(StringIO(bad))


def test_add_ions_requires_enough_water():
    mol = Xponge.load_pdb(StringIO("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"))
    with pytest.raises(ValueError, match="not enough WAT"):
        Xponge.Add_Ions(mol, {"NA": 1})
