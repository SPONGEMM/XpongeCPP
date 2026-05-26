from io import StringIO

import XpongeCPP as Xponge


def test_mol2_bond_order_refreshes_derived_state_for_gaff_assignment():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401

    assignment = Xponge.get_assignment_from_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
BEN
 12 12 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    C  -1.2131  -0.6884   0.0000   C.ar         1      BEN   0.000000
     2   C1  -1.2028   0.7064   0.0001   C.ar         1      BEN   0.000000
     3   C2  -0.0103  -1.3948   0.0000   C.ar         1      BEN   0.000000
     4   C3   0.0104   1.3948  -0.0001   C.ar         1      BEN   0.000000
     5   C4   1.2028  -0.7063   0.0000   C.ar         1      BEN   0.000000
     6   C5   1.2131   0.6884   0.0000   C.ar         1      BEN   0.000000
     7    H  -2.1577  -1.2244   0.0000   H            1      BEN   0.000000
     8   H1  -2.1393   1.2564   0.0001   H            1      BEN   0.000000
     9   H2  -0.0184  -2.4809  -0.0001   H            1      BEN   0.000000
    10   H3   0.0184   2.4808   0.0000   H            1      BEN   0.000000
    11   H4   2.1394  -1.2563   0.0001   H            1      BEN   0.000000
    12   H5   2.1577   1.2245   0.0000   H            1      BEN   0.000000
@<TRIPOS>BOND
     1      1      2 ar
     2      1      3 ar
     3      1      7 1
     4      2      4 ar
     5      2      8 1
     6      3      5 ar
     7      3      9 1
     8      4      6 ar
     9      4     10 1
    10      5      6 ar
    11      5     11 1
    12      6     12 1
@<TRIPOS>SUBSTRUCTURE
    1      BEN      1 ****               0 ****  ****
"""
        )
    )

    for atom in range(6):
        assert assignment.has_atom_marker(atom, "RG6")
        assert assignment.has_atom_marker(atom, "AR1")
    for atom1, atom2 in ((0, 1), (0, 2), (1, 3), (2, 4), (3, 5), (4, 5)):
        assert assignment.has_bond_marker(atom1, atom2, "AB")
        assert assignment.has_bond_marker(atom2, atom1, "AB")

    assignment.determine_atom_type("gaff")

    assert assignment.atom_types[:6] == ["ca"] * 6
    assert assignment.atom_types[6:] == ["ha"] * 6
