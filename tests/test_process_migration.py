import math
from io import StringIO

import numpy as np

import XpongeCPP as Xponge


def _load_chain():
    return Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
CHAIN
4 3 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
2 A2 1.0 0.0 0.0 C 1 LIG 0.0
3 A3 2.0 1.0 0.0 C 1 LIG 0.0
4 A4 3.0 1.0 1.0 C 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
"""
        )
    )


def _distance(atom1, atom2):
    return math.sqrt((atom1.x - atom2.x) ** 2 + (atom1.y - atom2.y) ** 2 + (atom1.z - atom2.z) ** 2)


def _angle(atom1, atom2, atom3):
    v1 = np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z])
    v2 = np.array([atom3.x - atom2.x, atom3.y - atom2.y, atom3.z - atom2.z])
    cosine = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
    cosine = max(-1.0, min(1.0, float(cosine)))
    return math.acos(cosine)


def _dihedral(atom1, atom2, atom3, atom4):
    p0 = np.array([atom1.x, atom1.y, atom1.z], dtype=float)
    p1 = np.array([atom2.x, atom2.y, atom2.z], dtype=float)
    p2 = np.array([atom3.x, atom3.y, atom3.z], dtype=float)
    p3 = np.array([atom4.x, atom4.y, atom4.z], dtype=float)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.atan2(y, x)


def test_impose_bond_updates_target_distance():
    mol = _load_chain()
    atom1, atom2 = mol.residues[0].atoms[:2]

    Xponge.impose_bond(mol, atom1, atom2, 2.0)

    assert math.isclose(_distance(atom1, atom2), 2.0, rel_tol=0.0, abs_tol=1e-6)


def test_impose_angle_updates_target_angle():
    mol = _load_chain()
    atom1, atom2, atom3 = mol.residues[0].atoms[:3]

    Xponge.impose_angle(mol, atom1, atom2, atom3, math.pi / 2)

    assert math.isclose(_angle(atom1, atom2, atom3), math.pi / 2, rel_tol=0.0, abs_tol=1e-6)


def test_impose_dihedral_updates_target_dihedral():
    mol = _load_chain()
    atom1, atom2, atom3, atom4 = mol.residues[0].atoms

    Xponge.impose_dihedral(mol, atom1, atom2, atom3, atom4, math.pi / 3)

    assert math.isclose(_dihedral(atom1, atom2, atom3, atom4), math.pi / 3, rel_tol=0.0, abs_tol=1e-6)


def test_main_axis_rotate_aligns_long_axis_to_requested_direction():
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
LINE
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
2 A2 5.0 0.0 0.0 C 1 LIG 0.0
3 A3 10.0 0.0 0.0 C 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""
        )
    )

    Xponge.main_axis_rotate(mol, direction_long=[0, 0, 1])

    xs = [atom.x for atom in mol.residues[0].atoms]
    zs = [atom.z for atom in mol.residues[0].atoms]
    assert max(xs) - min(xs) < 1e-6
    assert max(zs) - min(zs) > 9.9


def test_get_peptide_from_sequence_uses_terminal_templates():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.get_peptide_from_sequence("AG")

    assert [res.name for res in mol.residues] == ["NALA", "CGLY"]


def test_h_mass_repartition_moves_mass_from_heavy_atom_to_hydrogen():
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
CH
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C 1 LIG 0.0
2 H1 1.0 0.0 0.0 H 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )

    atom_c, atom_h = mol.residues[0].atoms
    Xponge.h_mass_repartition(mol, repartition_mass=1.1, repartition_rate=3)

    assert math.isclose(atom_h.mass, 3.024, rel_tol=0.0, abs_tol=1e-6)
    assert math.isclose(atom_c.mass, 9.994, rel_tol=0.0, abs_tol=1e-6)


def test_lattice_create_populates_block_region_with_simple_cubic_basis():
    basis = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
ATOM
1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
"""
        )
    )

    box = Xponge.BlockRegion(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, boundary=True)
    lattice = Xponge.Lattice("sc", basis_molecule=basis, scale=1.0)

    mol = lattice.create(box, box)

    assert mol.atom_count == 8
    assert mol.residue_count == 8
    assert list(mol.box_length) == [2.0, 2.0, 2.0]
