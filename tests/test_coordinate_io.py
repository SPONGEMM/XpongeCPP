from io import StringIO

import pytest
import XpongeCPP as Xponge


PDB_TEXT = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
END
"""


def test_load_coordinate_updates_molecule_and_box():
    mol = Xponge.load_pdb(StringIO(PDB_TEXT))
    text = """\
2
 1.0 2.0 3.0
 4.0 5.0 6.0
 10.0 11.0 12.0 80.0 90.0 100.0
"""

    coordinates, box = Xponge.load_coordinate(StringIO(text), mol)

    assert coordinates == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    assert box == [10.0, 11.0, 12.0, 80.0, 90.0, 100.0]
    assert mol.residues[0].atoms[0].x == 1.0
    assert mol.residues[0].atoms[1].z == 6.0
    assert mol.box_length == [10.0, 11.0, 12.0]
    assert mol.box_angle == [80.0, 90.0, 100.0]


def test_load_rst7_updates_molecule_and_box():
    mol = Xponge.load_pdb(StringIO(PDB_TEXT))
    text = """\
default name
    2
   1.0000000   2.0000000   3.0000000   4.0000000   5.0000000   6.0000000
  20.0000000  21.0000000  22.0000000  70.0000000  80.0000000  90.0000000
"""

    coordinates, box = Xponge.load_rst7(StringIO(text), mol)

    assert coordinates == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    assert box == [20.0, 21.0, 22.0, 70.0, 80.0, 90.0]
    assert mol.residues[0].atoms[0].y == 2.0
    assert mol.residues[0].atoms[1].x == 4.0
    assert mol.box_length == [20.0, 21.0, 22.0]
    assert mol.box_angle == [70.0, 80.0, 90.0]


def test_load_coordinate_accepts_three_value_box():
    text = """\
1
 1.0 2.0 3.0
 10.0 11.0 12.0
"""

    coordinates, box = Xponge.load_coordinate(StringIO(text))

    assert coordinates == [[1.0, 2.0, 3.0]]
    assert box == [10.0, 11.0, 12.0, 90.0, 90.0, 90.0]


def test_load_rst7_without_box_does_not_set_fake_box():
    mol = Xponge.load_pdb(StringIO(PDB_TEXT))
    text = """\
default name
    2
   1.0000000   2.0000000   3.0000000   4.0000000   5.0000000   6.0000000
"""

    coordinates, box = Xponge.load_rst7(StringIO(text), mol)

    assert coordinates == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    assert box == []
    assert mol.box_length == [0.0, 0.0, 0.0]
    assert mol.box_angle == [90.0, 90.0, 90.0]


def test_load_coordinate_rejects_atom_count_mismatch():
    mol = Xponge.load_pdb(StringIO(PDB_TEXT))
    text = """\
1
 1.0 2.0 3.0
 10.0 11.0 12.0
"""

    with pytest.raises(ValueError):
        Xponge.load_coordinate(StringIO(text), mol)
