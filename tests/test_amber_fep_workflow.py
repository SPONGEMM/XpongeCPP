from io import StringIO

import pytest
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


def test_merge_dual_topology_builds_a_b_ghost_atoms_and_match_map():
    import XpongeCPP.forcefield.special.fep as fep

    mol_a = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
A
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 ETH -0.2
2 H1 1.0 0.0 0.0 hc 1 ETH 0.1
3 H2 0.0 1.0 0.0 hc 1 ETH 0.1
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
"""
        )
    )
    mol_b = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
B
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 EOH -0.1
2 H1 1.0 0.0 0.0 hc 1 EOH 0.2
3 O1 0.0 1.0 0.0 oh 1 EOH -0.1
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
"""
        )
    )

    merged_from, merged_to, match_map = fep.merge_dual_topology(mol_a, 0, mol_b, {0: 0, 1: 1})

    assert match_map == {0: 0, 1: 1}
    assert merged_from.atom_count == 4
    assert merged_to.atom_count == 4
    assert merged_from.residues[0].name == "ETH_EOH"
    assert merged_to.residues[0].name == "EOH_ETH"
    from_atoms = merged_from.residues[0].atoms
    to_atoms = merged_to.residues[0].atoms
    assert [atom.name for atom in from_atoms] == ["C1", "H1", "H2", "O1R2"]
    assert from_atoms[3].type == "ZERO_LJ_ATOM"
    assert from_atoms[3].charge == 0.0
    assert from_atoms[3].subsys == 2
    assert to_atoms[2].type == "ZERO_LJ_ATOM"
    assert to_atoms[2].charge == 0.0
    assert to_atoms[2].subsys == 1
    assert to_atoms[3].type == "oh"


def test_merge_force_field_interpolates_charge_and_writes_softcore(tmp_path):
    import XpongeCPP.forcefield.special.fep as fep

    for atom_type, epsilon, rmin in [
        ("c3", 0.1094, 1.9080),
        ("hc", 0.0157, 1.4870),
        ("oh", 0.2100, 1.7210),
        ("ZERO_LJ_ATOM", 0.0, 0.0),
    ]:
        Xponge.register_amber_lj_parameter(atom_type, atom_type, epsilon, rmin)

    mol_a = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
A
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 ETH -0.2
2 H1 1.0 0.0 0.0 hc 1 ETH 0.2
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    mol_b = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
B
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 EOH 0.2
2 O1 1.0 0.0 0.0 oh 1 EOH -0.2
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    merged_from, merged_to, _ = fep.merge_dual_topology(mol_a, 0, mol_b, {0: 0})
    merged = fep.merge_force_field(merged_from, merged_to, 0.25, {"charge": 0.5})

    atoms = merged.residues[0].atoms
    assert atoms[0].charge == 0.0
    assert atoms[0].lj_type_b == "c3"
    assert atoms[1].charge == 0.1
    assert atoms[2].charge == -0.1
    assert atoms[1].subsys == 1
    assert atoms[2].subsys == 2

    out = Xponge.Save_SPONGE_Input(merged, prefix="merged", dirname=str(tmp_path))

    assert {"LJ_soft_core", "subsys_division", "bond_soft"}.issubset(out)
    assert "LJ" not in out
    assert (tmp_path / "merged_subsys_division.txt").read_text().splitlines() == ["3", "0", "1", "2"]
    assert (tmp_path / "merged_bond.txt").read_text().splitlines() == [
        "2",
        "0 1 225.000000 1.000000",
        "0 2 75.000000 1.000000",
    ]


def test_fep_compatibility_helpers_free_selected_residue_and_expose_xponge_names(tmp_path):
    import XpongeCPP.forcefield.special.fep as fep

    Xponge.register_amber_lj_parameter("c3", "c3", 0.1094, 1.9080)
    Xponge.register_amber_lj_parameter("ZERO_LJ_ATOM", "ZERO_LJ_ATOM", 0.0, 0.0)
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
FREE
1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 LIG -0.3
"""
        )
    )

    free = fep.get_free_molecule(mol, [0])
    assert free.residues[0].atoms[0].type == "ZERO_LJ_ATOM"
    assert free.residues[0].atoms[0].charge == 0.0
    assert free.residues[0].atoms[0].subsys == 1

    fep.save_soft_core_lj(free)
    out = Xponge.Save_SPONGE_Input(free, prefix="free", dirname=str(tmp_path))
    assert {"LJ_soft_core", "subsys_division"}.issubset(out)
    assert fep.Merge_Dual_Topology is fep.merge_dual_topology
    assert fep.Merge_Force_Field is fep.merge_force_field
    assert fep.Get_Free_Molecule is fep.get_free_molecule


def test_merge_dual_topology_can_derive_match_map_from_assign_when_rdkit_is_available():
    pytest.importorskip("rdkit")
    import XpongeCPP.forcefield.special.fep as fep

    mol_a = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
A
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 ETH 0.0
2 H1 1.0 0.0 0.0 hc 1 ETH 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    mol_b = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
B
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 EOH 0.0
2 H1 1.0 0.0 0.0 hc 1 EOH 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    assign_a = Xponge.get_assignment_from_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
A
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 ETH 0.0
2 H1 1.0 0.0 0.0 H 1 ETH 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    assign_b = Xponge.get_assignment_from_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
B
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 EOH 0.0
2 H1 1.0 0.0 0.0 H 1 EOH 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )

    merged_from, merged_to, match_map = fep.merge_dual_topology(mol_a, mol_a.residues[0], mol_b, assign_a, assign_b)

    assert match_map == {0: 0, 1: 1}
    assert merged_from.atom_count == merged_to.atom_count == 2
