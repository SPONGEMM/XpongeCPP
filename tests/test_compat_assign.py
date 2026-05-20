import XpongeCPP as Xponge


def test_legacy_atomtype_registry_supports_old_style_new_from_string_and_get_type():
    import Xponge

    Xponge.AtomType.New_From_String(
        """
name  mass   charge[e]
H     1.008  +0.25
C     12.00  -1.00
"""
    )

    h = Xponge.AtomType.get_type("H")
    c = Xponge.AtomType.Get_Type("C")
    all_types = Xponge.AtomType.get_all_types()

    assert h.name == "H"
    assert h.mass == 1.008
    assert h.__dict__["charge[e]"] == 0.25
    assert c.name == "C"
    assert c.mass == 12.0
    assert c.__dict__["charge[e]"] == -1.0
    assert all_types["H"] is h
    assert all_types["C"] is c


def test_legacy_basic_manual_building_shape_matches_old_test0_base_usage():
    import Xponge

    assign = Xponge.Assign()
    assign.add_atom("O", 0, 0, 0)
    assign.addAtom("H1", 0, 1, 0)
    assign.Add_Atom("H2", 1, 0, 0)
    assign.add_bond(0, 1, 1)
    assign.addBond(0, 2, 1)
    assign.Add_Bond(1, 2, 1)

    Xponge.AtomType.New_From_String(
        """
name  mass   charge[e]
H     1.008  +0.25
C     12.00  -1.00
"""
    )
    h = Xponge.AtomType.get_type("H")
    c = Xponge.AtomType.get_type("C")
    xyj = Xponge.ResidueType(name="XYJ")
    xyj.add_atom("H1", h, 1, 0, 0)
    xyj.add_atom("H2", "H", 0, 1, 0)
    xyj.add_atom("H3", h, -1, 0, 0)
    xyj.add_atom("H4", h, 0, -1, 0)
    xyj.add_atom("C", c, 0, 0, 0)
    mol = Xponge.Molecule(name="XYJ2")
    mol.add_residue(xyj)

    assert assign.atom_count == 3
    assert assign.bonds[0][1] == 1
    assert assign.bonds[1][2] == 1
    assert xyj.name == "XYJ"
    assert len(xyj.atoms) == 5
    assert mol.name == "XYJ2"
    assert mol.residue_count == 1


def test_legacy_assignrule_aliases_and_pure_string_path_match_old_usage():
    import Xponge
    from Xponge.assign import AssignRule

    rule = AssignRule("myrule", pure_string=True)

    @rule.Add_Rule("A", -1)
    def _a(i, assign):  # pylint: disable=unused-argument
        return True

    @rule.addRule("B", 1)
    def _b(i, assign):  # pylint: disable=unused-argument
        return True

    @rule.Set_Pre_Action
    def _pre(assign):
        assign.atoms[1] = "O"

    @rule.setPostAction
    def _post(assign):
        return assign

    assign = Xponge.Assign()
    assign.add_atom("H", 0, 0, 0)
    assign.add_atom("C", 1, 0, 0)
    assign.add_atom("O", 1, 0, 0)
    assign.add_bond(0, 1, 1)
    assign.add_bond(2, 1, 2)

    results = assign.determine_atom_type("myrule")
    assert results == ["B", "B", "B"]


def test_legacy_assign_delete_aliases_match_old_usage():
    import Xponge
    from io import StringIO

    s = StringIO(
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
    ben = Xponge.get_assignment_from_xyz(s)
    ben.add_bond(0, 11, 1)
    ben.add_atom("O", 0, 0, 0)
    ben.add_bond(0, 12, 1)
    ben.Delete_Atom(12)

    assert ben.atom_count == 12
    ben.deleteBond(0, 11)
    assert ben.bond_count == 12
    assert 11 not in ben.bonds[0]
    assert 0 not in ben.bonds[11]


def test_legacy_assign_common_camelcase_aliases_exist():
    import Xponge

    assign = Xponge.Assign()
    names = [
        "Determine_Ring_And_Bond_Type",
        "Save_As_PDB",
        "Set_Charge",
        "Set_Charges",
        "Set_Coordinate",
        "Set_Formal_Charge",
        "Set_Atom_Type",
        "Add_Bond_Marker",
        "Has_Bond_Marker",
    ]
    for name in names:
        assert hasattr(assign, name), name


def test_assign_atom_type_accepts_atomtype_objects_and_property_setter():
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")

    assign.set_atom_type(0, hw)
    assign.Set_Atom_Type(1, hw)
    assign.atom_types = [hw, hw, "HW", "HW"]

    assert assign.atom_types == ["HW", "HW", "HW", "HW"]


def test_assign_to_residuetype_handles_duplicate_atom_names_like_legacy_xponge():
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")
    assign.atom_types = [hw, hw, "HW", "HW"]

    residue_type = assign.to_residuetype("TES")

    assert [atom.name for atom in residue_type.atoms] == ["O", "O1", "H", "H1"]


def test_assign_to_residuetype_avoids_collisions_with_existing_numbered_names():
    import Xponge

    assign = Xponge.Assign("TES")
    assign.add_atom("C", 0.0, 0.0, 0.0, name="CA")
    assign.add_atom("C", 1.0, 0.0, 0.0, name="CA")
    assign.add_atom("C", 2.0, 0.0, 0.0, name="CA1")
    assign.add_bond(0, 1, 1)
    assign.add_bond(1, 2, 1)
    assign.atom_types = ["c3", "c3", "c3"]

    residue_type = assign.to_residuetype("TES")

    assert [atom.name for atom in residue_type.atoms] == ["CA", "CA1", "CA11"]


def test_assign_legacy_atom_numbers_and_residuetype_aliases_match_old_mokda_usage():
    import Xponge

    assign = Xponge.Assign("TES")

    assert assign.atom_numbers == 0
    assert callable(assign.toResidueType)
    assert callable(assign.To_ResidueType)

    assign.Add_Atom("ca", 0.0, 0.0, 0.0, "C1", 0.0)

    assert assign.atom_numbers == 1

    residue_type = assign.toResidueType("TES")

    assert residue_type.name == "TES"
    assert len(residue_type.atoms) == 1
