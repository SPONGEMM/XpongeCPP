from io import StringIO
from pathlib import Path

import pytest
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


CUSTOM_MOL2_TEXT = """\
@<TRIPOS>MOLECULE
CUSTOM_SOLVENT
5 3 2
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 O1 0.0000 0.0000 0.0000 OW 1 FAR -0.500
2 H1 8.0000 0.0000 0.0000 HW 1 FAR 0.250
3 H2 0.0000 8.0000 0.0000 HW 1 FAR 0.250
4 O2 2.0000 2.0000 2.0000 OW 2 LIG -0.500
5 H3 2.0000 2.0000 3.1000 HW 2 LIG 0.500
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
3 4 5 1
"""


FAR_WAT_MOL2_TEXT = """\
@<TRIPOS>MOLECULE
FAR_WAT_SYSTEM
8 5 3
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 O1 0.0000 0.0000 0.0000 OW 1 FAR -0.500
2 H1 8.0000 0.0000 0.0000 HW 1 FAR 0.500
3 O2 10.0000 0.0000 0.0000 OW 2 WAT -0.834
4 H2 10.9572 0.0000 0.0000 HW 2 WAT 0.417
5 H3 9.7610 0.9270 0.0000 HW 2 WAT 0.417
6 O3 14.0000 0.0000 0.0000 OW 3 WAT -0.834
7 H4 14.9572 0.0000 0.0000 HW 3 WAT 0.417
8 H5 13.7610 0.9270 0.0000 HW 3 WAT 0.417
@<TRIPOS>BOND
1 1 2 1
2 3 4 1
3 3 5 1
4 6 7 1
5 6 8 1
"""


def _snapshot_export_dir(path: Path) -> dict[str, str]:
    return {item.name: item.read_text() for item in sorted(path.iterdir()) if item.is_file()}


def _make_missing_oxt_pdb():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    template = Xponge.get_template_molecule("CLYS")
    lines = []
    serial = 1
    for atom in template.residues[0].atoms:
        if atom.name == "OXT":
            continue
        lines.append(
            f"ATOM  {serial:5d} {atom.name:>4s} CLYS A   1    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00          {atom.element:>2s}"
        )
        serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def test_load_pdb_preserves_molecule_residue_atom_layers():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.load_pdb(StringIO(PDB_TEXT))

    assert mol.name == "PDB"
    assert mol.atom_count == 4
    assert mol.residue_count == 1
    assert mol.residues[0].name == "NALA"
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


def test_add_solvent_box_appends_template_waters_and_exports_implicit_box(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    solute = Xponge.load_pdb(StringIO(PDB_TEXT))
    water = Xponge.load_mol2(StringIO(MOL2_TEXT))

    Xponge.Add_Solvent_Box(solute, water, 4.0, tolerance=2.5, n_solvent=4)

    assert solute.residue_count == 5
    assert solute.atom_count == 16
    assert solute.box_length[0] == 0.0
    assert solute.validate()
    Xponge.Save_SPONGE_Input(solute, prefix="case", dirname=str(tmp_path))
    box = [float(value) for value in (tmp_path / "case_coordinate.txt").read_text().splitlines()[-1].split()]
    assert box[0] > 0


def test_add_solvent_box_accepts_xponge_distance_vectors():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    scalar = Xponge.load_pdb(StringIO(PDB_TEXT))
    vector3 = Xponge.load_pdb(StringIO(PDB_TEXT))
    vector6 = Xponge.load_pdb(StringIO(PDB_TEXT))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(scalar, water, 4.0, tolerance=2.5)
    Xponge.Add_Solvent_Box(vector3, water, [4.0, 4.0, 4.0], tolerance=2.5)
    Xponge.Add_Solvent_Box(vector6, water, [0.0, 0.0, 0.0, 4.0, 4.0, 4.0], tolerance=2.5)

    assert vector3.residue_counts()["WAT"] == scalar.residue_counts()["WAT"]
    assert vector6.residue_counts()["WAT"] < scalar.residue_counts()["WAT"]
    with pytest.raises(TypeError):
        Xponge.Add_Solvent_Box(Xponge.load_pdb(StringIO(PDB_TEXT)), water, [1.0, 2.0], tolerance=2.5)


def test_add_solvent_box_replicates_all_residues_of_custom_solvent_unit(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    solute = Xponge.load_pdb(StringIO(PDB_TEXT))
    solvent = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    Xponge.Add_Solvent_Box(solute, solvent, 4.0, tolerance=2.5, n_solvent=1)
    Xponge.Save_SPONGE_Input(solute, prefix="customsol", dirname=str(tmp_path))

    assert solute.validate()
    assert solute.residue_count == 3
    assert [res.name for res in solute.residues] == ["NALA", "FAR", "LIG"]
    assert (tmp_path / "customsol_bond.txt").read_text().splitlines()[0] == "6"


def test_load_mol2_preserves_declared_residues_atom_properties():
    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    assert mol.name == "CUSTOM_SOLVENT"
    assert mol.residue_count == 2
    assert mol.atom_count == 5
    assert [res.name for res in mol.residues] == ["FAR", "LIG"]
    assert mol.residues[0].name2atom("O1").type == "OW"
    assert mol.residues[0].name2atom("H1").charge == 0.250
    assert mol.residues[1].name2atom("O2").type == "OW"
    assert mol.validate()


def test_load_mol2_as_template_registers_legacy_residuetype_lookup():
    Xponge.load_mol2(StringIO(MOL2_TEXT), as_template=True)

    wat = Xponge.ResidueType.get_type("WAT")

    assert wat.name == "WAT"
    assert Xponge.has_template("WAT")


def test_template_atoms_follow_forcefield_mass_and_element_inference():
    import importlib
    import XpongeCPP.forcefield.amber.ff14sb as ff14sb
    import XpongeCPP.forcefield.amber.tip3p as tip3p

    importlib.reload(ff14sb)
    importlib.reload(tip3p)

    hid = Xponge.get_template_molecule("HID")
    nd1 = hid.residues[0].name2atom("ND1")
    assert nd1.type == "NA"
    assert nd1.element == "N"
    assert nd1.mass == pytest.approx(14.01, abs=1e-6)

    sodium = Xponge.get_template_molecule("NA")
    assert sodium.residues[0].atoms[0].type == "Na+"
    assert sodium.residues[0].atoms[0].element == "Na"
    assert sodium.residues[0].atoms[0].mass == pytest.approx(22.99, abs=1e-6)


def test_assignment_from_residuetype_uses_mass_based_elements():
    hid = Xponge.ResidueType("HID_TMP")
    hid.add_atom("ND1", "NA", 0.0, 0.0, 0.0, charge=-0.38, mass=14.01)
    assignment = Xponge.get_assignment_from_residuetype(hid)

    assert assignment.atoms == ["N"]


def test_legacy_add_residue_link_accepts_atom_objects_and_atom_index_proxy():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    combined = Xponge.load_pdb(StringIO(PDB_TEXT)) + Xponge.load_pdb(StringIO(PDB_TEXT))
    atom_a = combined.residues[0].name2atom("C")
    atom_b = combined.residues[1].name2atom("N")

    combined.add_residue_link(atom_a, atom_b)

    assert combined.residue_links[-1] == [int(atom_a.index), int(atom_b.index)]
    assert combined.atom_index[atom_a] == int(atom_a.index)
    assert combined.atom_index[atom_b] == int(atom_b.index)


def test_residue_links_explicit_api_and_legacy_assignment_stay_in_sync(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    combined = Xponge.load_pdb(StringIO(PDB_TEXT)) + Xponge.load_pdb(StringIO(PDB_TEXT))
    atom_a = combined.residues[0].name2atom("C")
    atom_b = combined.residues[1].name2atom("N")
    combined.add_residue_link(atom_a, atom_b)

    saved_links = combined.get_residue_links()
    assert saved_links[-1] == [int(atom_a.index), int(atom_b.index)]

    combined.clear_residue_links()
    assert combined.residue_links == []

    cleared = tmp_path / "cleared.pdb"
    Xponge.Save_PDB(combined, cleared)
    assert "CONECT" not in cleared.read_text()

    combined.residue_links = saved_links
    assert combined.residue_links == saved_links

    restored = tmp_path / "restored.pdb"
    Xponge.Save_PDB(combined, restored)
    restored_text = restored.read_text()
    assert "CONECT" in restored_text
    assert f"{int(atom_a.index) + 1:5d}{int(atom_b.index) + 1:5d}" in restored_text


def test_set_residue_links_accepts_atom_objects_and_indices():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    combined = Xponge.load_pdb(StringIO(PDB_TEXT)) + Xponge.load_pdb(StringIO(PDB_TEXT))
    atom_a = combined.residues[0].name2atom("C")
    atom_b = combined.residues[1].name2atom("N")

    combined.set_residue_links([(atom_a, atom_b)])
    assert combined.residue_links == [[int(atom_a.index), int(atom_b.index)]]

    combined.residue_links = [(int(atom_a.index), int(atom_b.index))]
    assert combined.residue_links == [[int(atom_a.index), int(atom_b.index)]]


def test_add_missing_atoms_restores_terminal_oxt_without_moving_existing_atoms():
    import XpongeCPP.forcefield.amber as amber
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    data_dir = Path(__file__).resolve().parent / "data" / "8ryk"
    amber.load_parameters_from_frcmod(data_dir / "frcmod" / "interactive.frcmod", prefix=False)
    Xponge.load_mol2(data_dir / "edit_struct" / "CCS_3.gaff.mol2", as_template=True)
    mol = Xponge.load_pdb(
        str(data_dir / "8RYK_pdbfixer_H_ed.pdb"),
        ignore_conect=False,
        read_cryst1=False,
        unterminal_residues=["D:1"],
    )
    residue = mol.residues[-1]
    before = {atom.name: (atom.x, atom.y, atom.z) for atom in residue.atoms}

    mol.add_missing_atoms()

    residue = mol.residues[-1]
    after = {atom.name: (atom.x, atom.y, atom.z) for atom in residue.atoms}
    assert "OXT" in after
    for atom_name, coords in before.items():
        assert after[atom_name] == coords
    assert residue.name2atom("OXT").bad_coordinate is True


def test_load_mol2_declared_bonds_drive_topology_even_when_far(tmp_path):
    Xponge.register_tip3p()
    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    Xponge.Save_SPONGE_Input(mol, prefix="custom", dirname=str(tmp_path))

    assert (tmp_path / "custom_bond.txt").read_text().splitlines()[0] == "3"


def test_save_sponge_input_writes_core_files(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

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


def test_save_sponge_input_is_byte_identical_across_repeated_exports(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.load_pdb(StringIO(PDB_TEXT))
    water = Xponge.get_template_molecule("WAT")
    Xponge.Add_Solvent_Box(mol, water, 4.0, tolerance=2.5, n_solvent=4, seed=20260509)

    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    Xponge.Save_SPONGE_Input(mol, prefix="case", dirname=str(first))
    Xponge.Save_SPONGE_Input(mol, prefix="case", dirname=str(second))

    assert _snapshot_export_dir(first) == _snapshot_export_dir(second)


def test_save_sponge_input_writes_xponge_extra_bonded_force_files(tmp_path):
    Xponge.register_tip3p()
    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    mol.add_virtual_atom2(0, 1, 2, 3, 0.25, 0.75)
    mol.add_improper_dihedral(1, 2, 3, 4, 5.5, 180.0)
    mol.add_cmap_type(2, [9.0, 9.0, 9.0, 9.0])
    cmap_type = mol.add_cmap_type(2, [0.1, 0.2, 0.3, 0.4])
    mol.add_cmap(4, 3, 2, 1, 0, cmap_type)
    mol.add_nb14_extra(4, 1, 1.25, 2.5, 0.75)

    out = Xponge.Save_SPONGE_Input(mol, prefix="extra", dirname=str(tmp_path))

    assert {"virtual_atom", "improper_dihedral", "cmap", "nb14_extra"}.issubset(out)
    assert (tmp_path / "extra_virtual_atom.txt").read_text().splitlines() == [
        "2 0 1 2 3 0.250000 0.750000",
    ]
    assert (tmp_path / "extra_improper_dihedral.txt").read_text().splitlines() == [
        "1",
        "3 1 2 4 5.500000 180.000000",
    ]
    assert (tmp_path / "extra_nb14_extra.txt").read_text().splitlines() == [
        "1",
        "1 4 1.250000e+00 2.500000e+00 7.500000e-01",
    ]
    assert (tmp_path / "extra_cmap.txt").read_text().splitlines() == [
        "1 1",
        "2 ",
        "0.100000 0.200000 ",
        "0.300000 0.400000 ",
        "",
        "4 3 2 1 0 0",
    ]


def test_extra_bonded_force_indices_are_remapped_when_molecules_are_copied(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    Xponge.register_tip3p()
    target = Xponge.load_pdb(StringIO(PDB_TEXT))
    source = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))
    source.add_virtual_atom2(0, 1, 2, 3, 0.25, 0.75)
    source.add_improper_dihedral(1, 2, 3, 4, 5.5, 180.0)
    source.add_nb14_extra(4, 1, 1.25, 2.5, 0.75)
    cmap_type = source.add_cmap_type(2, [0.1, 0.2, 0.3, 0.4])
    source.add_cmap(4, 3, 2, 1, 0, cmap_type)

    Xponge.Add_Molecule(target, source)
    Xponge.Save_SPONGE_Input(target, prefix="merged_extra", dirname=str(tmp_path))

    assert target.validate()
    assert (tmp_path / "merged_extra_virtual_atom.txt").read_text().splitlines() == [
        "2 4 5 6 7 0.250000 0.750000",
    ]
    assert (tmp_path / "merged_extra_improper_dihedral.txt").read_text().splitlines() == [
        "1",
        "7 5 6 8 5.500000 180.000000",
    ]
    assert (tmp_path / "merged_extra_nb14_extra.txt").read_text().splitlines() == [
        "1",
        "5 8 1.250000e+00 2.500000e+00 7.500000e-01",
    ]
    assert (tmp_path / "merged_extra_cmap.txt").read_text().splitlines()[-1] == "8 7 6 5 4 0"


def test_extra_bonded_force_entries_on_removed_solvent_are_dropped_during_ion_replacement(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))
    mol.add_virtual_atom2(0, 1, 2, 3, 0.25, 0.75)
    mol.add_improper_dihedral(1, 2, 3, 4, 5.5, 180.0)
    mol.add_nb14_extra(4, 1, 1.25, 2.5, 0.75)

    Xponge.Add_Ions(mol, {"NA": 1}, seed=20260509, solvent="FAR")
    out = Xponge.Save_SPONGE_Input(mol, prefix="ion_extra", dirname=str(tmp_path))

    assert mol.validate()
    assert "virtual_atom" not in out
    assert "improper_dihedral" not in out
    assert "nb14_extra" not in out


def test_save_sponge_input_writes_xponge_general_bonded_force_files(tmp_path):
    Xponge.register_tip3p()
    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    mol.add_urey_bradley(2, 1, 0, 1.1, 2.2, 3.3, 4.4)
    mol.add_ryckaert_bellemans(4, 3, 2, 1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
    mol.add_bond_soft(4, 1, 6.5, 7.5, 2)

    out = Xponge.Save_SPONGE_Input(mol, prefix="general", dirname=str(tmp_path))

    assert {"urey_bradley", "Ryckaert_Bellemans", "bond_soft", "listed_forces"}.issubset(out)
    assert (tmp_path / "general_urey_bradley.txt").read_text().splitlines() == [
        "1",
        "0 1 2 1.100000 2.200000 3.300000 4.400000",
    ]
    assert (tmp_path / "general_Ryckaert_Bellemans.txt").read_text().splitlines() == [
        "1",
        "1 2 3 4 0.500000 1.500000 2.500000 3.500000 4.500000 5.500000",
    ]
    assert (tmp_path / "general_bond_soft.txt").read_text().splitlines() == [
        "1",
        "1 4 6.500000 7.500000 2",
    ]
    assert (tmp_path / "general_listed_forces.txt").read_text().startswith("[[[ Ryckaert_Bellemans ]]]")


def test_save_sponge_input_special_force_exports_are_byte_identical(tmp_path):
    Xponge.register_tip3p()
    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))
    mol.add_virtual_atom2(0, 1, 2, 3, 0.25, 0.75)
    mol.add_improper_dihedral(1, 2, 3, 4, 5.5, 180.0)
    mol.add_nb14_extra(4, 1, 1.25, 2.5, 0.75)
    cmap_type = mol.add_cmap_type(2, [0.1, 0.2, 0.3, 0.4])
    mol.add_cmap(4, 3, 2, 1, 0, cmap_type)
    mol.add_urey_bradley(2, 1, 0, 1.1, 2.2, 3.3, 4.4)
    mol.add_ryckaert_bellemans(4, 3, 2, 1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
    mol.add_bond_soft(4, 1, 6.5, 7.5, 2)

    first = tmp_path / "first_special"
    second = tmp_path / "second_special"
    first.mkdir()
    second.mkdir()
    Xponge.Save_SPONGE_Input(mol, prefix="special", dirname=str(first))
    Xponge.Save_SPONGE_Input(mol, prefix="special", dirname=str(second))

    assert _snapshot_export_dir(first) == _snapshot_export_dir(second)


def test_save_sponge_input_writes_xponge_special_state_files(tmp_path):
    Xponge.register_tip3p()
    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    mol.set_gb_radius("bondi_radii")
    mol.enable_min_bonded_parameters()
    mol.enable_subsys_division()
    mol.residues[0].atoms[0].subsys = 2
    mol.residues[0].atoms[1].bad_coordinate = True
    mol.residues[1].atoms[0].zero_lj_atom = True

    out = Xponge.Save_SPONGE_Input(mol, prefix="special", dirname=str(tmp_path))

    assert {"gb", "fake_mass", "fake_LJ", "fake_charge", "subsys_division"}.issubset(out)
    assert (tmp_path / "special_gb.txt").read_text().splitlines()[:3] == [
        "5",
        "1.5200 0.8500",
        "1.2000 0.8500",
    ]
    assert (tmp_path / "special_fake_mass.txt").read_text().splitlines() == [
        "5",
        "0.000",
        "1.000",
        "1.000",
        "1.000",
        "1.000",
    ]
    assert (tmp_path / "special_fake_charge.txt").read_text().splitlines() == [
        "5",
        "0.000000",
        "0.000000",
        "0.000000",
        "0.000000",
        "0.000000",
    ]
    assert (tmp_path / "special_subsys_division.txt").read_text().splitlines() == [
        "5",
        "2",
        "0",
        "0",
        "0",
        "0",
    ]
    assert (tmp_path / "special_fake_LJ.txt").read_text().splitlines()[0] == "5 1"


def test_set_gb_radius_matches_xponge_mass_based_rules():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    water = Xponge.load_mol2(StringIO(MOL2_TEXT))
    water.set_gb_radius("modified_bondi_radii")
    assert water.residues[0].name2atom("O").gb_radius == pytest.approx(1.5, abs=1e-6)
    assert water.residues[0].name2atom("H1").gb_radius == pytest.approx(0.8, abs=1e-6)
    assert water.residues[0].name2atom("H2").gb_radius == pytest.approx(0.8, abs=1e-6)

    peptide = Xponge.get_template_molecule("NALA")
    peptide.set_gb_radius("modified_bondi_radii")
    assert peptide.residues[0].name2atom("H1").gb_radius == pytest.approx(1.3, abs=1e-6)


def test_save_sponge_input_writes_xponge_pairwise_and_softcore_files(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))
    for residue in mol.residues:
        for atom in residue.atoms:
            atom.sw_type = "S"
            atom.edip_type = "E"
    mol.add_sw_type("S-S", 1.1, 2.2, 3.3, 4.0, 5.0, 6.6, 7.7, 8.8, 0.0, 0.0)
    mol.add_sw_type("S-S-S", 0.0, 0.0, 3.3, 0.0, 0.0, 0.0, 0.0, 0.0, 9.9, 10.1)
    mol.add_edip_type("E-E", 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 0.0, 0.0, 0.0, 12.12, 0.0, 0.0, 0.0, 0.0, 0.0)
    mol.add_edip_type("E-E-E", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.9, 10.1, 11.11, 12.12, 13.13, 0.0, 14.14, 15.15, 16.16, 17.17, 18.18)
    mol.enable_lj_soft_core()

    out = Xponge.Save_SPONGE_Input(mol, prefix="pairwise", dirname=str(tmp_path))

    assert {"SW", "EDIP", "LJ_soft_core", "subsys_division"}.issubset(out)
    assert "LJ" not in out
    assert (tmp_path / "pairwise_SW.txt").read_text().splitlines() == [
        "5 1",
        "# type1 type2 A B epsilon[kcal/mol] p q a gamma sigma[Angstrom] (This is the first required comment line)",
        "0 0 1.1 2.2 3.3 4.0 5.0 6.6 7.7 8.8",
        "# type1 type2 type3 lambda epsilon[kcal/mol] b (This is the second required comment line)",
        "0 0 0 9.9 3.3 10.1",
        "# atom type from the zeroth atom (This is the third required comment line)",
        "0",
        "0",
        "0",
        "0",
        "0",
    ]
    assert (tmp_path / "pairwise_EDIP.txt").read_text().splitlines()[:5] == [
        "5 1",
        "# type1 type2 alpha c[A] a[A] A[kcal/mol] B[A] rho beta sigma[A] (This is the first required comment line)",
        "0 0 5.5 4.4 3.3 1.1 2.2 0.0 6.6 12.12",
        "# type1 type2 type3 eta gamma[A] l[kcal/mol] Q0 mu u1 u2 u3 u4 (This is the second required comment line)",
        "0 0 0 9.9 10.1 11.11 14.14 12.12 15.15 16.16 17.17 18.18",
    ]
    assert (tmp_path / "pairwise_LJ_soft_core.txt").read_text().splitlines()[0] == "5 2 2"


def test_add_ions_randomly_replaces_waters_by_seed():
    Xponge.register_tip3p()
    solute = Xponge.load_pdb(StringIO(PDB_TEXT))
    water = Xponge.load_mol2(StringIO(MOL2_TEXT))
    Xponge.Add_Solvent_Box(solute, water, 4.0, tolerance=2.5, n_solvent=8)
    first_water_oxygen = tuple(
        getattr(solute.residues[1].atoms[0], axis) for axis in ("x", "y", "z")
    )

    Xponge.Add_Ions(solute, {"NA": 1}, seed=20260509)

    assert solute.residues[1].name == "NA"
    ion_position = tuple(getattr(solute.residues[1].atoms[0], axis) for axis in ("x", "y", "z"))
    assert ion_position != first_water_oxygen


def test_add_solvent_box_is_deterministic_for_fixed_seed():
    Xponge.register_tip3p()
    first = Xponge.load_pdb(StringIO(PDB_TEXT))
    second = Xponge.load_pdb(StringIO(PDB_TEXT))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(first, water, 4.0, tolerance=2.5, n_solvent=8, seed=20260509)
    Xponge.Add_Solvent_Box(second, water, 4.0, tolerance=2.5, n_solvent=8, seed=20260509)

    assert first.atom_count == second.atom_count
    assert first.residue_count == second.residue_count
    assert [res.name for res in first.residues] == [res.name for res in second.residues]
    assert [
        (atom.name, atom.type, atom.x, atom.y, atom.z)
        for atom in first.atoms
    ] == [
        (atom.name, atom.type, atom.x, atom.y, atom.z)
        for atom in second.atoms
    ]


def test_add_ions_preserves_mol2_explicit_bonds_after_rebuild(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.load_mol2(StringIO(FAR_WAT_MOL2_TEXT))

    Xponge.Add_Ions(mol, {"NA": 1}, seed=20260509)
    Xponge.Save_SPONGE_Input(mol, prefix="farwat", dirname=str(tmp_path))

    assert mol.validate()
    assert [res.name for res in mol.residues] == ["FAR", "NA", "WAT"]
    assert (tmp_path / "farwat_bond.txt").read_text().splitlines()[0] == "3"


def test_add_ions_can_replace_named_custom_solvent_residue(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    Xponge.Add_Ions(mol, {"NA": 1}, seed=20260509, solvent="FAR")
    Xponge.Save_SPONGE_Input(mol, prefix="customion", dirname=str(tmp_path))

    assert mol.validate()
    assert [res.name for res in mol.residues] == ["LIG", "NA"]
    assert (tmp_path / "customion_bond.txt").read_text().splitlines()[0] == "1"


def test_add_molecule_preserves_source_explicit_bonds(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    Xponge.register_tip3p()
    target = Xponge.load_pdb(StringIO(PDB_TEXT))
    source = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT))

    Xponge.Add_Molecule(target, source)
    Xponge.Save_SPONGE_Input(target, prefix="merged", dirname=str(tmp_path))

    assert target.validate()
    assert target.atom_count == 9
    assert target.residue_count == 3
    assert [res.name for res in target.residues] == ["NALA", "FAR", "LIG"]
    assert (tmp_path / "merged_bond.txt").read_text().splitlines()[0] == "6"


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
