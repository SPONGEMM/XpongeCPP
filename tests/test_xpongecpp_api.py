from io import StringIO

import XpongeCPP as Xponge


PDB_TEXT = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.300   2.360   0.000  1.00  0.00           O
END
"""


MOL2_TEXT = """\
@<TRIPOS>MOLECULE
WAT
3 2 1
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 O 0.0000 0.0000 0.0000 OW 1 WAT -0.834
2 H1 0.9572 0.0000 0.0000 HW 1 WAT 0.417
3 H2 -0.2390 0.9270 0.0000 HW 1 WAT 0.417
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
"""


def test_load_pdb_preserves_molecule_residue_atom_layers():
    mol = Xponge.load_pdb(StringIO(PDB_TEXT))

    assert mol.name == "PDB"
    assert mol.atom_count == 4
    assert mol.residue_count == 1
    assert mol.residues[0].name == "ALA"
    assert mol.residues[0].atom_count == 4
    assert mol.residues[0].name2atom("CA").name == "CA"


def test_residue_type_is_writable_and_versioned():
    restype = Xponge.ResidueType("TMP")
    assert restype.version == 0

    restype.add_atom("A", "C", 0.0, 0.0, 0.0, charge=0.1, mass=12.0)
    restype.add_atom("B", "H", 1.0, 0.0, 0.0, charge=-0.1, mass=1.0)
    restype.add_connectivity("A", "B")

    assert restype.version == 3
    assert restype.atom_count == 2
    assert restype.bond_count == 1


def test_add_solvent_box_appends_template_waters_and_sets_box():
    solute = Xponge.load_pdb(StringIO(PDB_TEXT))
    water = Xponge.load_mol2(StringIO(MOL2_TEXT))

    Xponge.Add_Solvent_Box(solute, water, 4.0, tolerance=2.5, n_solvent=4)

    assert solute.residue_count == 5
    assert solute.atom_count == 16
    assert solute.box_length[0] > 0
    assert solute.validate()


def test_save_sponge_input_writes_core_files(tmp_path):
    mol = Xponge.load_pdb(StringIO(PDB_TEXT))
    Xponge.Set_Box_Padding(mol, 3.0)

    out = Xponge.Save_SPONGE_Input(mol, prefix="case", dirname=str(tmp_path))

    assert sorted(out) == [
        "LJ",
        "angle",
        "atom_name",
        "atom_type_name",
        "bond",
        "charge",
        "coordinate",
        "dihedral",
        "exclude",
        "mass",
        "nb14",
        "residue",
        "resname",
    ]
    assert (tmp_path / "case_coordinate.txt").exists()
    assert (tmp_path / "case_residue.txt").read_text().splitlines()[0] == "4 1"


def test_assign_builds_graph_markers_and_residue_type():
    assign = Xponge.Assign("ASN")
    assign.add_atom("O", 0.0, 0.0, 0.0, name="O")
    assign.add_atom("H", 0.96, 0.0, 0.0, name="H1")
    assign.add_atom("H", -0.24, 0.93, 0.0, name="H2")
    assign.determine_connectivity(simple_cutoff=1.2)
    assign.determine_atom_type("amber")

    assert assign.atom_types == ["O", "H", "H"]
    restype = assign.to_residuetype("WATX")
    assert restype.atom_count == 3
    assert restype.bond_count == 2
