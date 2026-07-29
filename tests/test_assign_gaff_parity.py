from __future__ import annotations

import json
import subprocess
import sys
from base64 import b64encode
from io import StringIO
from pathlib import Path
from textwrap import dedent

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import XpongeCPP as Xponge
from XpongeCPP.helper.rdkit import rdmol_to_assign

ORIGIN_REPO = Path("/mnt/data8t/Software/Xponge/Xponge-origin")


def _assigned_type_names(assign, indices: list[int]) -> list[str]:
    return [str(assign.atom_types[i]) for i in indices]


def _assignment_from_smiles_with_3d(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260603
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        status = AllChem.EmbedMolecule(mol, randomSeed=20260603, useRandomCoords=True)
    assert status == 0
    if AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    return rdmol_to_assign(mol)


def _origin_reference_types_from_smiles(smiles: str, rule: str, indices: list[int]) -> list[str]:
    if not ORIGIN_REPO.exists():
        pytest.skip("local Xponge-origin repo not available")
    env_site_packages = Path(sys.executable).resolve().parent.parent / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages"
    script = dedent(
        f"""
        import json, sys
        sys.path.insert(0, {str(ORIGIN_REPO)!r})
        sys.path.append({str(env_site_packages)!r})
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import Xponge
        from Xponge.helper.rdkit import rdmol_to_assign
        import Xponge.forcefield.amber.gaff
        import Xponge.forcefield.amber.gaff2
        mol = Chem.MolFromSmiles({smiles!r})
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 20260603
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            status = AllChem.EmbedMolecule(mol, randomSeed=20260603, useRandomCoords=True)
        assert status == 0
        if AllChem.UFFHasAllMoleculeParams(mol):
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        assign = rdmol_to_assign(mol)
        assign.determine_atom_type({rule!r})
        print(json.dumps([assign.atom_types[i].name for i in {indices!r}]))
        """
    )
    output = subprocess.check_output([sys.executable, "-S", "-c", script], text=True)
    payload = output.strip().splitlines()[-1]
    return json.loads(payload)


def _origin_reference_types_from_mol2_text(mol2_text: str, rule: str, indices: list[int]) -> list[str]:
    if not ORIGIN_REPO.exists():
        pytest.skip("local Xponge-origin repo not available")
    env_site_packages = (
        Path(sys.executable).resolve().parent.parent
        / "lib"
        / f"python{sys.version_info.major}.{sys.version_info.minor}"
        / "site-packages"
    )
    encoded = b64encode(mol2_text.encode("utf-8")).decode("ascii")
    script = dedent(
        f"""
        import base64, json, sys
        from io import StringIO
        sys.path.insert(0, {str(ORIGIN_REPO)!r})
        sys.path.append({str(env_site_packages)!r})
        import Xponge
        import Xponge.forcefield.amber.gaff
        import Xponge.forcefield.amber.gaff2
        text = base64.b64decode({encoded!r}).decode("utf-8")
        assign = Xponge.get_assignment_from_mol2(StringIO(text))
        assign.determine_atom_type({rule!r})
        print(json.dumps([assign.atom_types[i].name for i in {indices!r}]))
        """
    )
    output = subprocess.check_output([sys.executable, "-S", "-c", script], text=True)
    payload = output.strip().splitlines()[-1]
    return json.loads(payload)


def test_gaff_cross_family_alternating_regression():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401

    assign = _assignment_from_smiles_with_3d("C#CC1Cc2ccsc2C(N)=N1")
    assign.determine_atom_type("gaff")
    assert _assigned_type_names(assign, [4, 5, 6, 8, 9]) == ["cc", "cc", "cd", "cd", "cf"]


def test_public_get_assignment_from_smiles_supports_gaff_and_gaff2():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("c1ccccc1")
    assign.determine_atom_type("gaff")
    assert _assigned_type_names(assign, list(range(6))) == ["ca"] * 6
    assert _assigned_type_names(assign, list(range(6, 12))) == ["ha"] * 6

    assign = Xponge.get_assignment_from_smiles("c1ccccc1")
    assign.determine_atom_type("gaff2")
    assert _assigned_type_names(assign, list(range(6))) == ["ca"] * 6
    assert _assigned_type_names(assign, list(range(6, 12))) == ["ha"] * 6


def test_gaff_atom_type_determination():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401

    s = StringIO(
        dedent(
            """12
            BEN
            C -1.213  -0.688   0.000
            C -1.203   0.706   0.000
            C -0.010  -1.395   0.000
            C  0.010   1.395  -0.000
            C  1.213   0.688   0.000
            C  1.203  -0.706   0.000
            H  0.018   2.481   0.000
            H -2.158  -1.224   0.000
            H -2.139   1.256   0.000
            H -0.018  -2.481  -0.000
            H  2.139  -1.256   0.000
            H  2.158   1.224   0.000
            """
        )
    )
    ben = Xponge.get_assignment_from_xyz(s)
    ben.determine_atom_type("gaff")
    assert _assigned_type_names(ben, list(range(6))) == ["ca"] * 6
    assert _assigned_type_names(ben, list(range(6, 12))) == ["ha"] * 6


def test_gaff2_atom_type_determination():
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    s = StringIO(
        dedent(
            """12
            BEN
            C -1.213  -0.688   0.000
            C -1.203   0.706   0.000
            C -0.010  -1.395   0.000
            C  0.010   1.395  -0.000
            C  1.213   0.688   0.000
            C  1.203  -0.706   0.000
            H  0.018   2.481   0.000
            H -2.158  -1.224   0.000
            H -2.139   1.256   0.000
            H -0.018  -2.481  -0.000
            H  2.139  -1.256   0.000
            H  2.158   1.224   0.000
            """
        )
    )
    ben = Xponge.get_assignment_from_xyz(s)
    ben.determine_atom_type("gaff2")
    assert _assigned_type_names(ben, list(range(6))) == ["ca"] * 6
    assert _assigned_type_names(ben, list(range(6, 12))) == ["ha"] * 6


def test_gaff_assignment_state_and_residuetype_regression():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    mol2_text = dedent(
        """
            @<TRIPOS>MOLECULE
            *****
             49 52 0 0 0
            SMALL
            GASTEIGER

            @<TRIPOS>ATOM
                  1 C          -5.0279    0.1289   -3.6674 C.3     1  UNL1        0.0790
                  2 O          -5.3788   -0.4045   -2.3945 O.3     1  UNL1       -0.4914
                  3 C          -4.5173   -0.4542   -1.2868 C.ar    1  UNL1        0.1635
                  4 C          -3.2000    0.0321   -1.3447 C.ar    1  UNL1        0.0099
                  5 C          -2.3695   -0.0330   -0.2234 C.ar    1  UNL1        0.0816
                  6 N          -1.1016    0.4431   -0.3066 N.ar    1  UNL1       -0.2141
                  7 C          -0.2561    0.4015    0.7552 C.ar    1  UNL1        0.2201
                  8 N           1.0848    0.9190    0.6348 N.pl3   1  UNL1       -0.2966
                  9 C           2.0519    0.7649    1.7317 C.3     1  UNL1        0.0307
                 10 C           3.4509    0.4867    1.1717 C.3     1  UNL1        0.0309
                 11 N           3.8109    1.4396    0.1093 N.am    1  UNL1       -0.2920
                 12 C           5.1991    1.6674   -0.1941 C.2     1  UNL1        0.2820
                 13 O           5.6311    2.8484   -0.2818 O.2     1  UNL1       -0.2666
                 14 C           6.1395    0.5544   -0.4037 C.ar    1  UNL1        0.1892
                 15 C           5.8198   -0.7390   -0.7717 C.ar    1  UNL1       -0.0113
                 16 C           7.0220   -1.4055   -0.8700 C.ar    1  UNL1       -0.0231
                 17 C           7.9992   -0.4810   -0.5636 C.ar    1  UNL1        0.0928
                 18 O           7.4427    0.6775   -0.2948 O.2     1  UNL1       -0.4583
                 19 C           2.8100    2.3820   -0.4034 C.3     1  UNL1        0.0309
                 20 C           1.4526    1.6893   -0.5640 C.3     1  UNL1        0.0307
                 21 N          -0.6894   -0.1269    1.9299 N.ar    1  UNL1       -0.1992
                 22 C          -1.9460   -0.6220    2.0847 C.ar    1  UNL1        0.1291
                 23 N          -2.3033   -1.1577    3.3656 N.pl3   1  UNL1       -0.3425
                 24 C          -2.8413   -0.5883    0.9865 C.ar    1  UNL1        0.0437
                 25 C          -4.1632   -1.0751    1.0392 C.ar    1  UNL1       -0.0038
                 26 C          -5.0002   -1.0089   -0.0932 C.ar    1  UNL1        0.1623
                 27 O          -6.3216   -1.4849   -0.0730 O.3     1  UNL1       -0.4914
                 28 C          -6.9550   -2.0733    1.0587 C.3     1  UNL1        0.0790
                 29 H          -4.7609    1.2025   -3.5712 H       1  UNL1        0.0660
                 30 H          -5.8976    0.0393   -4.3500 H       1  UNL1        0.0660
                 31 H          -4.1776   -0.4404   -4.0988 H       1  UNL1        0.0660
                 32 H          -2.8103    0.4633   -2.2575 H       1  UNL1        0.0677
                 33 H           1.7741   -0.0780    2.4004 H       1  UNL1        0.0480
                 34 H           2.0682    1.6966    2.3374 H       1  UNL1        0.0480
                 35 H           3.4693   -0.5455    0.7672 H       1  UNL1        0.0480
                 36 H           4.1850    0.5516    2.0045 H       1  UNL1        0.0480
                 37 H           4.8408   -1.1449   -0.9847 H       1  UNL1        0.0657
                 38 H           7.1705   -2.4403   -1.1478 H       1  UNL1        0.0649
                 39 H           9.0649   -0.6656   -0.5476 H       1  UNL1        0.1029
                 40 H           2.7137    3.2346    0.3036 H       1  UNL1        0.0480
                 41 H           3.1255    2.7820   -1.3922 H       1  UNL1        0.0480
                 42 H           1.4933    1.0028   -1.4374 H       1  UNL1        0.0480
                 43 H           0.6833    2.4672   -0.7626 H       1  UNL1        0.0480
                 44 H          -3.2383   -1.5553    3.5907 H       1  UNL1        0.1437
                 45 H          -1.5965   -1.1556    4.1344 H       1  UNL1        0.1437
                 46 H          -4.5564   -1.5055    1.9429 H       1  UNL1        0.0663
                 47 H          -6.9926   -1.3457    1.8968 H       1  UNL1        0.0660
                 48 H          -7.9929   -2.3547    0.7870 H       1  UNL1        0.0660
                 49 H          -6.4089   -2.9890    1.3691 H       1  UNL1        0.0660
            @<TRIPOS>BOND
                 1     1     2    1
                 2     2     3    1
                 3     3     4   ar
                 4     4     5   ar
                 5     5     6   ar
                 6     6     7   ar
                 7     7     8    1
                 8     8     9    1
                 9     9    10    1
                10    10    11    1
                11    11    12   am
                12    12    13    2
                13    12    14    1
                14    14    15   ar
                15    15    16   ar
                16    16    17   ar
                17    17    18   ar
                18    11    19    1
                19    19    20    1
                20     7    21   ar
                21    21    22   ar
                22    22    23    1
                23    22    24   ar
                24    24    25   ar
                25    25    26   ar
                26    26    27    1
                27    27    28    1
                28    26     3   ar
                29    24     5   ar
                30    20     8    1
                31    18    14   ar
                32     1    29    1
                33     1    30    1
                34     1    31    1
                35     4    32    1
                36     9    33    1
                37     9    34    1
                38    10    35    1
                39    10    36    1
                40    15    37    1
                41    16    38    1
                42    17    39    1
                43    19    40    1
                44    19    41    1
                45    20    42    1
                46    20    43    1
                47    23    44    1
                48    23    45    1
                49    25    46    1
                50    28    47    1
                51    28    48    1
                52    28    49    1
        """
    )
    s = StringIO(mol2_text)
    assign = Xponge.get_assignment_from_mol2(s)
    assert assign.built is False

    assign.determine_atom_type("gaff2")
    assert _assigned_type_names(assign, [13, 14, 15, 16]) == ["cc", "cd", "cd", "cc"]
    assert _assigned_type_names(assign, [13, 14, 15, 16]) == _origin_reference_types_from_mol2_text(
        mol2_text, "gaff2", [13, 14, 15, 16]
    )

    assign.determine_atom_type("gaff")
    assert _assigned_type_names(assign, [13, 14, 15, 16]) == ["cc", "cd", "cd", "cc"]
    assert _assigned_type_names(assign, [13, 14, 15, 16]) == _origin_reference_types_from_mol2_text(
        mol2_text, "gaff", [13, 14, 15, 16]
    )

    restype = assign.to_residuetype("PRZ")
    assert [str(atom.type) for atom in restype.atoms[13:17]] == ["cc", "cd", "cd", "cc"]


def test_gaff_cp_cq_pure_aromatic_regression():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    smiles = "c1ccc(-c2nc(-c3ccccc3)c(-c3ccccc3)nc2-c2ccccc2)cc1"
    expected = ["cp", "cp", "cp", "cp", "cp", "cp", "cq", "cq"]

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff")
    assert _assigned_type_names(assign, [3, 4, 6, 7, 13, 14, 21, 22]) == expected

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff2")
    assert _assigned_type_names(assign, [3, 4, 6, 7, 13, 14, 21, 22]) == expected


def test_gaff_nitroso_ne_regression():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    smiles = "CN(C)S(=O)(=O)c1cc2c(N=O)c(O)[nH]c2c2c1CCCC2"
    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff")
    assert str(assign.atom_types[10]) == "ne"

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff2")
    assert str(assign.atom_types[10]) == "ne"


def test_gaff_nitroso_sequence_sensitive_n2_regression():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    smiles = "O=Nc1c(O)[nH]c2c3c(c([N+](=O)[O-])cc12)CCCC3"
    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff")
    assert str(assign.atom_types[1]) == "n2"

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff2")
    assert str(assign.atom_types[1]) == "n2"


def test_gaff_carbonyl_ring_c_regression():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    smiles = "Cc1noc2c(-c3ccccc3)nn(CCCN3CCN(c4cccc(Cl)c4)CC3)c(=O)c12"
    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff")
    assert str(assign.atom_types[30]) == "c"
    assert str(assign.atom_types[32]) == "cc"

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff2")
    assert str(assign.atom_types[30]) == "c"
    assert str(assign.atom_types[32]) == "cc"


def test_gaff2_special_nitrogen_types():
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    s = StringIO(
        dedent(
            """
            @<TRIPOS>MOLECULE
            MOL
               14    13     1     0     0
            SMALL
            No Charge or Current Charge

            @<TRIPOS>ATOM
                  1 C1          3.5400    1.4200    0.0000 C.3       1 MOL     -0.444095
                  2 H1          2.5460    1.2490   -0.3830 H         1 MOL      0.214061
                  3 H2          4.2190    0.6730   -0.3830 H         1 MOL      0.214061
                  4 H3          3.5320    1.3950    1.0800 H         1 MOL      0.214061
                  5 N1          4.0060    2.7710   -0.4390 N.4       1 MOL      0.063003
                  6 H4          4.0060    2.7710   -1.4480 H         1 MOL      0.342733
                  7 C2          3.0680    3.8500    0.0010 C.3       1 MOL     -0.444095
                  8 H5          3.4180    4.7960   -0.3830 H         1 MOL      0.214061
                  9 H6          2.0830    3.6360   -0.3830 H         1 MOL      0.214061
                 10 H7          3.0510    3.8700    1.0800 H         1 MOL      0.214061
                 11 C3          5.4100    3.0440    0.0010 C.3       1 MOL     -0.444095
                 12 H8          5.7170    4.0040   -0.3830 H         1 MOL      0.214061
                 13 H9          6.0550    2.2680   -0.3830 H         1 MOL      0.214061
                 14 H10         5.4360    3.0480    1.0800 H         1 MOL      0.214061
            @<TRIPOS>BOND
                 1    2    1 1
                 2    3    1 1
                 3    4    1 1
                 4    5    1 1
                 5    6    5 1
                 6    7    5 1
                 7    8    7 1
                 8    9    7 1
                 9   10    7 1
                10   11    5 1
                11   12   11 1
                12   13   11 1
                13   14   11 1
            @<TRIPOS>SUBSTRUCTURE
                 1 MOL         1 TEMP              0 ****  ****    0 ROOT
            """
        )
    )
    me3nh = Xponge.get_assignment_from_mol2(s)
    me3nh.determine_atom_type("gaff2")
    assert str(me3nh.atom_types[4]) == "nx"
    assert str(me3nh.atom_types[5]) == "hn"


def test_gaff2_amide_and_thiocarbonyl_types():
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    nma = StringIO(
        dedent(
            """
            @<TRIPOS>MOLECULE
            NMA
                9     8     1     0     0
            SMALL
            USER_CHARGES
            @<TRIPOS>ATOM
                  1 C1          0.0000    0.0000    0.0000 C.3       1 NMA      0.000000
                  2 H1         -0.3600   -0.5800    0.8900 H         1 NMA      0.000000
                  3 H2         -0.3600   -0.5800   -0.8900 H         1 NMA      0.000000
                  4 H3         -0.5400    0.9800    0.0000 H         1 NMA      0.000000
                  5 C2          1.5200   -0.0400    0.0000 C.2       1 NMA      0.000000
                  6 O1          2.1800   -1.0700    0.0000 O.2       1 NMA      0.000000
                  7 N1          2.2500    1.0900    0.0000 N.am      1 NMA      0.000000
                  8 H4          1.7600    1.9600    0.0000 H         1 NMA      0.000000
                  9 C3          3.6900    1.1800    0.0000 C.3       1 NMA      0.000000
            @<TRIPOS>BOND
                 1     1     2 1
                 2     1     3 1
                 3     1     4 1
                 4     1     5 1
                 5     5     6 2
                 6     5     7 1
                 7     7     8 1
                 8     7     9 1
            @<TRIPOS>SUBSTRUCTURE
                 1 NMA         1 TEMP              0 ****  ****    0 ROOT
            """
        )
    )
    nma = Xponge.get_assignment_from_mol2(nma)
    nma.determine_atom_type("gaff2")
    assert str(nma.atom_types[6]) == "ns"

    thiocarbonyl = StringIO(
        dedent(
            """
            @<TRIPOS>MOLECULE
            THI
                4     3     1     0     0
            SMALL
            USER_CHARGES
            @<TRIPOS>ATOM
                  1 C1          0.0000    0.0000    0.0000 C.2       1 THI      0.000000
                  2 S1          1.6200    0.0000    0.0000 S.2       1 THI      0.000000
                  3 H1         -0.5800    0.9300    0.0000 H         1 THI      0.000000
                  4 H2         -0.5800   -0.9300    0.0000 H         1 THI      0.000000
            @<TRIPOS>BOND
                 1     1     2 2
                 2     1     3 1
                 3     1     4 1
            @<TRIPOS>SUBSTRUCTURE
                 1 THI         1 TEMP              0 ****  ****    0 ROOT
            """
        )
    )
    thiocarbonyl = Xponge.get_assignment_from_mol2(thiocarbonyl)
    thiocarbonyl.determine_atom_type("gaff2")
    assert str(thiocarbonyl.atom_types[0]) == "cs"


def test_gaff2_ring_sp3_carbons():
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    cyclopentane = StringIO(
        dedent(
            """
            @<TRIPOS>MOLECULE
            CYP
               15    15     1     0     0
            SMALL
            USER_CHARGES
            @<TRIPOS>ATOM
                  1 C1          1.2140    0.0000    0.0000 C.3       1 CYP      0.000000
                  2 C2          0.3750    1.1540    0.0000 C.3       1 CYP      0.000000
                  3 C3         -0.9820    0.7130    0.0000 C.3       1 CYP      0.000000
                  4 C4         -0.9820   -0.7130    0.0000 C.3       1 CYP      0.000000
                  5 C5          0.3750   -1.1540    0.0000 C.3       1 CYP      0.000000
                  6 H1          2.2940    0.0000    0.0000 H         1 CYP      0.000000
                  7 H2          0.7880    1.9970    0.0000 H         1 CYP      0.000000
                  8 H3         -1.8510    1.3440    0.0000 H         1 CYP      0.000000
                  9 H4         -1.8510   -1.3440    0.0000 H         1 CYP      0.000000
                 10 H5          0.7880   -1.9970    0.0000 H         1 CYP      0.000000
                 11 H6          1.2140    0.0000    1.0900 H         1 CYP      0.000000
                 12 H7          0.3750    1.1540    1.0900 H         1 CYP      0.000000
                 13 H8         -0.9820    0.7130    1.0900 H         1 CYP      0.000000
                 14 H9         -0.9820   -0.7130    1.0900 H         1 CYP      0.000000
                 15 H10         0.3750   -1.1540    1.0900 H         1 CYP      0.000000
            @<TRIPOS>BOND
                 1     1     2 1
                 2     2     3 1
                 3     3     4 1
                 4     4     5 1
                 5     5     1 1
                 6     1     6 1
                 7     2     7 1
                 8     3     8 1
                 9     4     9 1
                10     5    10 1
                11     1    11 1
                12     2    12 1
                13     3    13 1
                14     4    14 1
                15     5    15 1
            @<TRIPOS>SUBSTRUCTURE
                 1 CYP         1 TEMP              0 ****  ****    0 ROOT
            """
        )
    )
    cyclopentane = Xponge.get_assignment_from_mol2(cyclopentane)
    cyclopentane.determine_atom_type("gaff2")
    assert _assigned_type_names(cyclopentane, list(range(5))) == ["c5"] * 5

    cyclohexane = StringIO(
        dedent(
            """
            @<TRIPOS>MOLECULE
            CYH
               18    18     1     0     0
            SMALL
            USER_CHARGES
            @<TRIPOS>ATOM
                  1 C1          1.2140    0.7010    0.0000 C.3       1 CYH      0.000000
                  2 C2          0.0000    1.4020    0.0000 C.3       1 CYH      0.000000
                  3 C3         -1.2140    0.7010    0.0000 C.3       1 CYH      0.000000
                  4 C4         -1.2140   -0.7010    0.0000 C.3       1 CYH      0.000000
                  5 C5          0.0000   -1.4020    0.0000 C.3       1 CYH      0.000000
                  6 C6          1.2140   -0.7010    0.0000 C.3       1 CYH      0.000000
                  7 H1          2.1570    1.2450    0.0000 H         1 CYH      0.000000
                  8 H2          1.2140    0.7010    1.0900 H         1 CYH      0.000000
                  9 H3          0.0000    2.4900    0.0000 H         1 CYH      0.000000
                 10 H4          0.0000    1.4020    1.0900 H         1 CYH      0.000000
                 11 H5         -2.1570    1.2450    0.0000 H         1 CYH      0.000000
                 12 H6         -1.2140    0.7010    1.0900 H         1 CYH      0.000000
                 13 H7         -2.1570   -1.2450    0.0000 H         1 CYH      0.000000
                 14 H8         -1.2140   -0.7010    1.0900 H         1 CYH      0.000000
                 15 H9          0.0000   -2.4900    0.0000 H         1 CYH      0.000000
                 16 H10         0.0000   -1.4020    1.0900 H         1 CYH      0.000000
                 17 H11         2.1570   -1.2450    0.0000 H         1 CYH      0.000000
                 18 H12         1.2140   -0.7010    1.0900 H         1 CYH      0.000000
            @<TRIPOS>BOND
                 1     1     2 1
                 2     2     3 1
                 3     3     4 1
                 4     4     5 1
                 5     5     6 1
                 6     6     1 1
                 7     1     7 1
                 8     1     8 1
                 9     2     9 1
                10     2    10 1
                11     3    11 1
                12     3    12 1
                13     4    13 1
                14     4    14 1
                15     5    15 1
                16     5    16 1
                17     6    17 1
                18     6    18 1
            @<TRIPOS>SUBSTRUCTURE
                 1 CYH         1 TEMP              0 ****  ****    0 ROOT
            """
        )
    )
    cyclohexane = Xponge.get_assignment_from_mol2(cyclohexane)
    cyclohexane.determine_atom_type("gaff2")
    assert _assigned_type_names(cyclohexane, list(range(6))) == ["c6"] * 6


def test_gaff2_sulfoxide_s4_regression():
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    assign = _assignment_from_smiles_with_3d("C[S+]([O-])c1ccccc1")
    assign.determine_atom_type("gaff2")
    assert str(assign.atom_types[1]) == "s4"
    assert str(assign.atom_types[3]) == "ca"


def test_gaff_macrocycle_matches_origin_reference():
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    smiles = (
        "C=CC1=C(C)c2cc3nc(cc4[nH]c(cc5[nH]c(cc1n2)c(C)c5CCC(=O)NC1C(O)OC(CO)C(O)C1O)"
        "c(CCC(=O)NC1C(O)OC(CO)C(O)C1O)c4C)C(C=C)=C3C"
    )
    indices = [6, 8, 10, 14, 18, 20, 21, 23, 40, 57]

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff")
    assert _assigned_type_names(assign, indices) == _origin_reference_types_from_smiles(smiles, "gaff", indices)

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type("gaff2")
    assert _assigned_type_names(assign, indices) == _origin_reference_types_from_smiles(smiles, "gaff2", indices)


@pytest.mark.parametrize(
    ("rule", "smiles", "indices"),
    [
        ("gaff", "C#CC1Cc2ccsc2C(N)=N1", [4, 5, 6, 8, 9]),
        ("gaff", "c1ccc(-c2nc(-c3ccccc3)c(-c3ccccc3)nc2-c2ccccc2)cc1", [3, 4, 6, 7, 13, 14, 21, 22]),
        ("gaff2", "c1ccc(-c2nc(-c3ccccc3)c(-c3ccccc3)nc2-c2ccccc2)cc1", [3, 4, 6, 7, 13, 14, 21, 22]),
        ("gaff", "CN(C)S(=O)(=O)c1cc2c(N=O)c(O)[nH]c2c2c1CCCC2", [10]),
        ("gaff2", "CN(C)S(=O)(=O)c1cc2c(N=O)c(O)[nH]c2c2c1CCCC2", [10]),
        ("gaff", "O=Nc1c(O)[nH]c2c3c(c([N+](=O)[O-])cc12)CCCC3", [1]),
        ("gaff2", "O=Nc1c(O)[nH]c2c3c(c([N+](=O)[O-])cc12)CCCC3", [1]),
        ("gaff", "Cc1noc2c(-c3ccccc3)nn(CCCN3CCN(c4cccc(Cl)c4)CC3)c(=O)c12", [30, 32]),
        ("gaff2", "Cc1noc2c(-c3ccccc3)nn(CCCN3CCN(c4cccc(Cl)c4)CC3)c(=O)c12", [30, 32]),
        ("gaff2", "C[S+]([O-])c1ccccc1", [1, 3]),
    ],
)
def test_curated_smiles_cases_match_local_origin_reference(rule, smiles, indices):
    module = "gaff2" if rule == "gaff2" else "gaff"
    __import__(f"XpongeCPP.forcefield.amber.{module}")

    assign = _assignment_from_smiles_with_3d(smiles)
    assign.determine_atom_type(rule)
    assert _assigned_type_names(assign, indices) == _origin_reference_types_from_smiles(smiles, rule, indices)


@pytest.mark.parametrize(
    ("rule", "mol2_text", "indices"),
    [
        (
            "gaff2",
            dedent(
                """
                @<TRIPOS>MOLECULE
                MOL
                   14    13     1     0     0
                SMALL
                No Charge or Current Charge

                @<TRIPOS>ATOM
                      1 C1          3.5400    1.4200    0.0000 C.3       1 MOL     -0.444095
                      2 H1          2.5460    1.2490   -0.3830 H         1 MOL      0.214061
                      3 H2          4.2190    0.6730   -0.3830 H         1 MOL      0.214061
                      4 H3          3.5320    1.3950    1.0800 H         1 MOL      0.214061
                      5 N1          4.0060    2.7710   -0.4390 N.4       1 MOL      0.063003
                      6 H4          4.0060    2.7710   -1.4480 H         1 MOL      0.342733
                      7 C2          3.0680    3.8500    0.0010 C.3       1 MOL     -0.444095
                      8 H5          3.4180    4.7960   -0.3830 H         1 MOL      0.214061
                      9 H6          2.0830    3.6360   -0.3830 H         1 MOL      0.214061
                     10 H7          3.0510    3.8700    1.0800 H         1 MOL      0.214061
                     11 C3          5.4100    3.0440    0.0010 C.3       1 MOL     -0.444095
                     12 H8          5.7170    4.0040   -0.3830 H         1 MOL      0.214061
                     13 H9          6.0550    2.2680   -0.3830 H         1 MOL      0.214061
                     14 H10         5.4360    3.0480    1.0800 H         1 MOL      0.214061
                @<TRIPOS>BOND
                     1    2    1 1
                     2    3    1 1
                     3    4    1 1
                     4    5    1 1
                     5    6    5 1
                     6    7    5 1
                     7    8    7 1
                     8    9    7 1
                     9   10    7 1
                    10   11    5 1
                    11   12   11 1
                    12   13   11 1
                    13   14   11 1
                @<TRIPOS>SUBSTRUCTURE
                     1 MOL         1 TEMP              0 ****  ****    0 ROOT
                """
            ),
            [4, 5],
        ),
        (
            "gaff2",
            dedent(
                """
                @<TRIPOS>MOLECULE
                THI
                    4     3     1     0     0
                SMALL
                USER_CHARGES
                @<TRIPOS>ATOM
                      1 C1          0.0000    0.0000    0.0000 C.2       1 THI      0.000000
                      2 S1          1.6200    0.0000    0.0000 S.2       1 THI      0.000000
                      3 H1         -0.5800    0.9300    0.0000 H         1 THI      0.000000
                      4 H2         -0.5800   -0.9300    0.0000 H         1 THI      0.000000
                @<TRIPOS>BOND
                     1     1     2 2
                     2     1     3 1
                     3     1     4 1
                @<TRIPOS>SUBSTRUCTURE
                     1 THI         1 TEMP              0 ****  ****    0 ROOT
                """
            ),
            [0],
        ),
    ],
)
def test_curated_mol2_cases_match_local_origin_reference(rule, mol2_text, indices):
    module = "gaff2" if rule == "gaff2" else "gaff"
    __import__(f"XpongeCPP.forcefield.amber.{module}")

    assign = Xponge.get_assignment_from_mol2(StringIO(mol2_text))
    assign.determine_atom_type(rule)
    assert _assigned_type_names(assign, indices) == _origin_reference_types_from_mol2_text(mol2_text, rule, indices)
