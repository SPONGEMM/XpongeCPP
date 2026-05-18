import pytest


def test_addsolventbox_camelcase_alias_executes_real_process_workflow(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    monkeypatch.chdir(tmp_path)

    mol = Xponge.ACE + Xponge.ALA + Xponge.NME
    solvated = Xponge.AddSolventBox(mol, Xponge.WAT, 8, n_solvent=2)
    outputs = Xponge.Save_SPONGE_Input(solvated, "solv")

    assert outputs is solvated
    assert solvated.residue_count >= 2
    assert solvated.residues[0].name == "ACE"
    assert (tmp_path / "solv_coordinate.txt").is_file()
    assert (tmp_path / "solv_resname.txt").is_file()


def test_process_camelcase_aliases_execute_real_geometry_and_mass_workflow():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    import math

    bond_mol = Xponge.ACE + Xponge.ALA + Xponge.NME
    bond_before = math.dist(
        (bond_mol.atoms[0].x, bond_mol.atoms[0].y, bond_mol.atoms[0].z),
        (bond_mol.atoms[1].x, bond_mol.atoms[1].y, bond_mol.atoms[1].z),
    )
    Xponge.Impose_Bond(bond_mol, bond_mol.atoms[0], bond_mol.atoms[1], 1.5)
    bond_after = math.dist(
        (bond_mol.atoms[0].x, bond_mol.atoms[0].y, bond_mol.atoms[0].z),
        (bond_mol.atoms[1].x, bond_mol.atoms[1].y, bond_mol.atoms[1].z),
    )
    assert bond_after != bond_before

    assign = Xponge.Assign()
    for i, (x, y, z) in enumerate(((0, 0, 0), (1, 0, 0), (2, 1, 0), (3, 1, 1)), start=1):
        assign.add_atom(f"A{i}", x, y, z)
    assign.add_bond(0, 1, 1)
    assign.add_bond(1, 2, 1)
    assign.add_bond(2, 3, 1)
    for i in range(4):
        assign.set_atom_type(i, "C")
    residue_type = assign.to_residuetype("LIN")

    angle_mol = Xponge.Molecule(name="LINA")
    angle_mol.add_residue(residue_type)
    angle_before = (float(angle_mol.atoms[2].x), float(angle_mol.atoms[2].y), float(angle_mol.atoms[2].z))
    Xponge.Impose_Angle(angle_mol, angle_mol.atoms[0], angle_mol.atoms[1], angle_mol.atoms[2], 2.1)
    angle_after = (float(angle_mol.atoms[2].x), float(angle_mol.atoms[2].y), float(angle_mol.atoms[2].z))
    assert angle_after != angle_before

    dihedral_mol = Xponge.Molecule(name="LIND")
    dihedral_mol.add_residue(residue_type)
    dihedral_before = (float(dihedral_mol.atoms[3].x), float(dihedral_mol.atoms[3].y), float(dihedral_mol.atoms[3].z))
    Xponge.Impose_Dihedral(
        dihedral_mol,
        dihedral_mol.atoms[0],
        dihedral_mol.atoms[1],
        dihedral_mol.atoms[2],
        dihedral_mol.atoms[3],
        1.1,
    )
    dihedral_after = (float(dihedral_mol.atoms[3].x), float(dihedral_mol.atoms[3].y), float(dihedral_mol.atoms[3].z))
    assert dihedral_after != dihedral_before

    mass_mol = Xponge.ACE + Xponge.ALA + Xponge.NME
    masses_before = [float(atom.mass) for atom in mass_mol.atoms]
    Xponge.H_Mass_Repartition(mass_mol)
    masses_after = [float(atom.mass) for atom in mass_mol.atoms]
    assert masses_after != masses_before
    assert pytest.approx(sum(masses_after), rel=1e-12, abs=1e-12) == sum(masses_before)


def test_process_sequence_box_and_solvent_replace_aliases_execute_real_workflow(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    peptide = Xponge.Get_Peptide_From_Sequence("AGA")
    assert [res.name for res in peptide.residues] == ["NALA", "GLY", "CALA"]

    before_box = list(peptide.box_length)
    Xponge.SetBoxPadding(peptide, 5.0)
    after_box = list(peptide.box_length)
    assert after_box != before_box

    solvated = Xponge.AddSolventBox(peptide, Xponge.WAT, 8, n_solvent=3)
    before_names = [res.name for res in solvated.residues]
    Xponge.Solvent_Replace(solvated, Xponge.WAT, {Xponge.NA: 1}, sort=True)
    after_names = [res.name for res in solvated.residues]

    assert before_names.count("NA") == 0
    assert after_names.count("NA") == 1
    assert after_names.count("WAT") == before_names.count("WAT") - 1

    monkeypatch.chdir(tmp_path)
    outputs = Xponge.Save_SPONGE_Input(solvated, "replace")
    assert outputs is solvated
    assert (tmp_path / "replace_resname.txt").is_file()


def test_main_axis_rotate_alias_executes_real_process_workflow():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.ACE + Xponge.ALA + Xponge.NME
    before = (
        float(mol.atoms[0].x),
        float(mol.atoms[0].y),
        float(mol.atoms[0].z),
    )
    Xponge.Main_Axis_Rotate(mol)
    after = (
        float(mol.atoms[0].x),
        float(mol.atoms[0].y),
        float(mol.atoms[0].z),
    )
    assert after != before


def test_add_molecule_alias_executes_real_process_workflow():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.ACE
    other = Xponge.ALA + Xponge.NME

    Xponge.Add_Molecule(mol, other)

    assert mol.residue_count == 3
    assert [res.name for res in mol.residues] == ["ACE", "ALA", "NME"]


def test_add_ions_alias_accepts_template_like_keys_and_executes():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.AddSolventBox(Xponge.ACE + Xponge.ALA + Xponge.NME, Xponge.WAT, 8, n_solvent=3)
    before_names = [res.name for res in mol.residues]

    Xponge.AddIons(mol, {Xponge.NA: 1})

    after_names = [res.name for res in mol.residues]
    assert after_names.count("NA") == before_names.count("NA") + 1
