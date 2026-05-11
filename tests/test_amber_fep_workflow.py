from io import StringIO

import XpongeCPP as Xponge


def test_fep_module_prepares_softcore_state_and_exports_expected_files(tmp_path):
    import XpongeCPP.forcefield.special.fep as fep

    Xponge.register_amber_lj_parameter("c3", "c3", 0.1094, 1.9080)
    Xponge.register_amber_lj_parameter("hc", "hc", 0.0157, 1.4870)
    Xponge.register_amber_lj_parameter("oh", "oh", 0.2100, 1.7210)

    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
FEP
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 LIG 0.0
2 H1 1.0 0.0 0.0 hc 1 LIG 0.0
3 O1 2.0 0.0 0.0 oh 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
"""
        )
    )

    fep.prepare_lj_soft_core(mol, {0: "c3", 1: "hc", 2: "oh"}, subsys=1)
    out = Xponge.Save_SPONGE_Input(mol, prefix="fep", dirname=str(tmp_path))

    assert {"LJ_soft_core", "subsys_division"}.issubset(out)
    assert "LJ" not in out
    assert (tmp_path / "fep_subsys_division.txt").read_text().splitlines() == ["3", "1", "1", "1"]


def test_fep_bonded_merge_helpers_follow_xponge_soft_bond_shape(tmp_path):
    import XpongeCPP.forcefield.special.fep as fep

    Xponge.register_amber_lj_parameter("A", "A", 0.1, 1.0)
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
FEPBOND
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 A 1 LIG 0.0
2 A2 1.0 0.0 0.0 A 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    fep.add_soft_bond_from_a(mol, 0, 1, k=10.0, b=1.5)
    fep.add_soft_bond_from_b(mol, 0, 1, k=20.0, b=1.8)

    out = Xponge.Save_SPONGE_Input(mol, prefix="soft", dirname=str(tmp_path))

    assert "bond_soft" in out
    assert (tmp_path / "soft_bond_soft.txt").read_text().splitlines() == [
        "2",
        "0 1 10.000000 1.500000 0",
        "0 1 20.000000 1.800000 1",
    ]
