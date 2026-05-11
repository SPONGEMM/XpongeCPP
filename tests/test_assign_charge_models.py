import math
import io

import pytest
import XpongeCPP as Xponge


def _assignment(name, atoms, bonds, formal_charges=None):
    assignment = Xponge.Assign(name)
    for index, element in enumerate(atoms):
        assignment.add_atom(element, float(index), 0.0, 0.0, f"{element}{index + 1}", 0.0)
    for atom1, atom2, order in bonds:
        assignment.add_bond(atom1, atom2, order)
    if formal_charges:
        for atom, charge in formal_charges.items():
            assignment.set_formal_charge(atom, charge)
    return assignment


def _assert_charges_close(actual, expected, tol=1e-6):
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        assert math.isclose(got, want, abs_tol=tol)


def test_tpacm4_calculate_charge_matches_xponge_methane_and_ethane():
    methane = _assignment(
        "methane",
        ["C", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
    )
    methane.calculate_charge("tpacm4")
    _assert_charges_close(methane.charges, [-0.28116, 0.07029, 0.07029, 0.07029, 0.07029])
    assert math.isclose(sum(methane.charges), 0.0, abs_tol=1e-9)

    ethane = _assignment(
        "ethane",
        ["C", "C", "H", "H", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)],
    )
    ethane.calculate_charge("TPACM4")
    _assert_charges_close(
        ethane.charges,
        [-0.205298, -0.205298, 0.068433, 0.068433, 0.068433, 0.068433, 0.068433, 0.068433],
    )
    assert math.isclose(sum(ethane.charges), 0.0, abs_tol=1e-9)


def test_tpacm4_calculate_charge_uses_aromatic_and_functional_group_context():
    benzene = _assignment(
        "benzene",
        ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 2),
            (4, 5, 1),
            (5, 0, 2),
            (0, 6, 1),
            (1, 7, 1),
            (2, 8, 1),
            (3, 9, 1),
            (4, 10, 1),
            (5, 11, 1),
        ],
    )
    benzene.calculate_charge("tpacm4")
    _assert_charges_close(benzene.charges[:6], [-0.155438] * 6)
    _assert_charges_close(benzene.charges[6:], [0.155438] * 6)

    methyl_acetate = _assignment(
        "methyl_acetate",
        ["C", "C", "O", "O", "C", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (0, 5, 1),
            (0, 6, 1),
            (0, 7, 1),
            (1, 2, 2),
            (1, 3, 1),
            (3, 4, 1),
            (4, 8, 1),
            (4, 9, 1),
            (4, 10, 1),
        ],
    )
    methyl_acetate.calculate_charge("tpacm4")
    _assert_charges_close(
        methyl_acetate.charges,
        [-0.403173, 0.812999, -0.551684, -0.366511, -0.018544, 0.086684, 0.086684, 0.086684, 0.088952, 0.088952, 0.088952],
    )


def test_tpacm4_respects_explicit_total_charge_and_formal_charge_property():
    methyl_ammonium = _assignment(
        "methyl_ammonium",
        ["C", "N", "H", "H", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)],
        {1: 1},
    )
    assert methyl_ammonium.formal_charges == [0, 1, 0, 0, 0, 0, 0, 0]

    methyl_ammonium.calculate_charge("tpacm4")
    _assert_charges_close(
        methyl_ammonium.charges,
        [-0.046932, -0.223509, 0.116603, 0.116603, 0.116603, 0.306877, 0.306877, 0.306877],
    )
    assert math.isclose(sum(methyl_ammonium.charges), 1.0, abs_tol=1e-9)


def test_calculate_charge_reports_optional_dependencies_clearly(monkeypatch):
    import builtins

    assignment = _assignment("methane", ["C", "H", "H", "H", "H"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    original_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("rdkit") or name.startswith("pyscf"):
            raise ImportError(name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="RDKit"):
        assignment.calculate_charge("gasteiger")
    with pytest.raises(ImportError, match="PySCF"):
        assignment.calculate_charge("resp")


def test_assign_charge_aliases_and_pubchem_signature_are_xponge_compatible(monkeypatch):
    assignment = _assignment("methane", ["C", "H", "H", "H", "H"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assignment.Calculate_Charge("tpacm4")
    _assert_charges_close(assignment.charges, [-0.28116, 0.07029, 0.07029, 0.07029, 0.07029])

    import builtins

    original_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("pubchempy"):
            raise ImportError(name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="PubChemPy"):
        Xponge.get_assignment_from_pubchem("CC", "smiles")


def test_cif_assignment_returns_lattice_info_like_xponge():
    cif = """
data_demo
_cell_length_a    10.0
_cell_length_b    11.0
_cell_length_c    12.0
_cell_angle_alpha 90.0
_cell_angle_beta  91.0
_cell_angle_gamma 92.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_Cartn_x
_atom_site_Cartn_y
_atom_site_Cartn_z
C1 C 0.0 0.0 0.0
H1 H 1.0 0.0 0.0
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif)
    assert assignment.name == "demo"
    assert assignment.atoms == ["C", "H"]
    assert lattice["cell_length"] == [10.0, 11.0, 12.0]
    assert lattice["cell_angle"] == [90.0, 91.0, 92.0]


def test_set_ph_deprotonates_carboxylic_acid_like_xponge_phmodel():
    acetic_acid = _assignment(
        "acetic_acid",
        ["C", "C", "O", "O", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 4, 1), (0, 5, 1), (0, 6, 1), (1, 2, 2), (1, 3, 1), (3, 7, 1)],
    )

    total_charge = acetic_acid.set_ph(7.0)

    assert total_charge == -1
    assert acetic_acid.atoms == ["C", "C", "O", "O", "H", "H", "H"]
    assert acetic_acid.formal_charges == [0, 0, 0, -1, 0, 0, 0]
    assert acetic_acid.bond_count == 6


def test_determine_connectivity_accepts_xponge_tolerance_signature():
    assignment = Xponge.Assign("water")
    assignment.add_atom("O", 0.0, 0.0, 0.0)
    assignment.add_atom("H", 0.96, 0.0, 0.0)
    assignment.add_atom("H", -0.24, 0.93, 0.0)

    assignment.determine_connectivity(tolerance=1.0)

    assert assignment.bond_count == 2
    assert sorted((min(i, j), max(i, j), order) for i, bonds in enumerate(assignment.bonds) for j, order in bonds.items() if i < j) == [
        (0, 1, -1),
        (0, 2, -1),
    ]


def test_pdb_and_xyz_assignment_entrypoints_accept_xponge_parameters():
    pdb = """\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000                       C
HETATM    2  H1  LIG A   1       1.090   0.000   0.000                       H
HETATM    3  C2  SOL A   2       5.000   0.000   0.000                       C
HETATM    4  H2  SOL A   2       6.090   0.000   0.000                       H
END
"""
    assignment = Xponge.get_assignment_from_pdb(io.StringIO(pdb), only_residue="LIG", bond_tolerance=1.0, total_charge=0)
    assert assignment.name == "LIG"
    assert assignment.atoms == ["C", "H"]
    assert assignment.bond_count == 1

    xyz = """2
fragment
C 0.000 0.000 0.000
H 1.090 0.000 0.000
"""
    assignment = Xponge.get_assignment_from_xyz(io.StringIO(xyz), bond_tolerance=1.0, total_charge=0)
    assert assignment.name == "fragment"
    assert assignment.bond_count == 1


def test_cif_assignment_supports_fractional_coordinates_and_symmetry_basis():
    cif = """
data_frac
_cell_length_a    10.0
_cell_length_b    20.0
_cell_length_c    30.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
_symmetry_equiv_pos_as_xyz
'x,y,z
x+1/2,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.5 0.25 0.1
H1 H 0.6 0.25 0.1
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 H1 1.0
"""
    assignment, lattice = Xponge.get_assignment_from_cif(cif, keep_cell_angle=False)
    assert assignment.name == "frac"
    assert assignment.coordinates[0] == pytest.approx([5.0, 5.0, 3.0])
    assert assignment.bond_count == 1
    assert lattice["cell_angle"] == [90, 90, 90]
    assert lattice["basis_position"] == [[1, 1, 1], [1.5, 1, 1]]


def test_phmodel_atom_type_and_set_ph_cover_phenol_and_alcohol():
    phenol = _assignment(
        "phenol",
        ["C", "C", "C", "C", "C", "C", "O", "H", "H", "H", "H", "H", "H"],
        [
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 2),
            (4, 5, 1),
            (5, 0, 2),
            (0, 7, 1),
            (1, 8, 1),
            (2, 9, 1),
            (3, 10, 1),
            (4, 11, 1),
            (5, 6, 1),
            (6, 12, 1),
        ],
    )
    types = phenol.determine_atom_type("phmodel")
    assert types[12] == "A-phenol"
    assert phenol.set_ph(11.0) == -1
    assert phenol.atoms.count("H") == 5

    alcohol = _assignment("methanol", ["C", "O", "H", "H", "H", "H"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1)])
    assert alcohol.determine_atom_type("phmodel")[5] == "A-alcohol"


def test_save_as_mol2_atomtype_argument_and_equal_atoms_api(tmp_path):
    ethane = _assignment(
        "ethane",
        ["C", "C", "H", "H", "H", "H", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)],
    )
    ethane.calculate_charge("tpacm4")
    mol2 = tmp_path / "ethane.mol2"
    ethane.save_as_mol2(str(mol2), atomtype="sybyl")
    text = mol2.read_text()
    assert "@<TRIPOS>UNITY_ATOM_ATTR" not in text

    with pytest.raises(ValueError, match="atomtype"):
        ethane.save_as_mol2(str(tmp_path / "bad.mol2"), atomtype="gaff")

    try:
        groups = ethane.determine_equal_atoms()
    except ImportError as exc:
        pytest.skip(str(exc))
    assert any(set(group) == {0, 1} for group in groups)


def test_resp_uses_pyscf_backend_or_reports_missing_dependency():
    water = _assignment("water", ["O", "H", "H"], [(0, 1, 1), (0, 2, 1)])
    try:
        water.calculate_charge("resp", basis="sto-3g", charge=0, grid_density=1, grid_cell_layer=1, only_esp=True)
    except ImportError as exc:
        assert "PySCF" in str(exc)
        return
    assert len(water.charges) == 3
    assert math.isclose(sum(water.charges), 0.0, abs_tol=1e-5)
