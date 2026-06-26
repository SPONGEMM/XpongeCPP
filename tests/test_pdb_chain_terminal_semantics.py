from io import StringIO

import pytest

import XpongeCPP as Xponge


def _load_amber():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401


def _bond_pairs(path):
    pairs = set()
    for line in path.read_text().splitlines()[1:]:
        if line.strip():
            words = line.split()
            pairs.add(tuple(sorted((int(words[0]), int(words[1])))))
    return pairs


def _residue_atom(residue, atom_name):
    name2atom = (
        getattr(residue, "name2atom", None)
        or getattr(residue, "Name2Atom", None)
        or getattr(residue, "Name2atom", None)
    )
    return name2atom(atom_name) if callable(name2atom) else getattr(residue, atom_name)


def test_pdb_ter_breaks_topology_chain_link(tmp_path):
    _load_amber()
    pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       7.000   0.000   0.000  1.00  0.00           O
TER
ATOM      9  N   ALA B   1       0.000   5.000   0.000  1.00  0.00           N
ATOM     10  CA  ALA B   1       1.000   5.000   0.000  1.00  0.00           C
ATOM     11  C   ALA B   1       2.000   5.000   0.000  1.00  0.00           C
ATOM     12  O   ALA B   1       3.000   5.000   0.000  1.00  0.00           O
ATOM     13  N   GLY B   2       4.000   5.000   0.000  1.00  0.00           N
ATOM     14  CA  GLY B   2       5.000   5.000   0.000  1.00  0.00           C
ATOM     15  C   GLY B   2       6.000   5.000   0.000  1.00  0.00           C
ATOM     16  O   GLY B   2       7.000   5.000   0.000  1.00  0.00           O
TER
END
"""
    mol = Xponge.load_pdb(StringIO(pdb))

    Xponge.Save_SPONGE_Input(mol, prefix="chain", dirname=str(tmp_path))
    pairs = _bond_pairs(tmp_path / "chain_bond.txt")

    assert [res.name for res in mol.residues] == ["NALA", "CGLY", "NALA", "CGLY"]
    assert (2, 4) in pairs
    assert (10, 12) in pairs
    assert (6, 8) not in pairs


def test_pdb_unterminal_residues_match_xponge_selectors():
    _load_amber()
    pdb = """\
ATOM      1  N   ALA A  12       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A  12       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A  12       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A  12       3.000   0.000   0.000  1.00  0.00           O
TER
ATOM      5  N   GLY B  35A      0.000   5.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY B  35A      1.000   5.000   0.000  1.00  0.00           C
ATOM      7  C   GLY B  35A      2.000   5.000   0.000  1.00  0.00           C
ATOM      8  O   GLY B  35A      3.000   5.000   0.000  1.00  0.00           O
TER
END
"""

    by_int = Xponge.load_pdb(StringIO(pdb), unterminal_residues=[12])
    by_string = Xponge.load_pdb(StringIO(pdb), unterminal_residues=["B:35A"])
    by_tuple = Xponge.load_pdb(StringIO(pdb), unterminal_residues=[("B", 35, "A")])

    assert [res.name for res in by_int.residues] == ["ALA", "NGLY"]
    assert [res.name for res in by_string.residues] == ["NALA", "GLY"]
    assert [res.name for res in by_tuple.residues] == ["NALA", "GLY"]


def test_pdb_explicit_terminal_residues_disable_inference():
    _load_amber()
    pdb = """\
ATOM      1  N   VAL A   2       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  VAL A   2       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   VAL A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   VAL A   2       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   TRP A   3       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  TRP A   3       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   TRP A   3       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   TRP A   3       7.000   0.000   0.000  1.00  0.00           O
TER
"""
    mol = Xponge.load_pdb(
        StringIO(pdb),
        ignore_hydrogen=True,
        terminal_residues=[{"chain_id": "A", "res_seq": 2, "n_terminal": True}],
        infer_terminals=False,
    )

    assert [res.name for res in mol.residues] == ["NVAL", "TRP"]

    blank_chain_pdb = """\
ATOM      1  N   VAL     2       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  VAL     2       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   VAL     2       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   VAL     2       3.000   0.000   0.000  1.00  0.00           O
TER
"""
    blank_chain = Xponge.load_pdb(
        StringIO(blank_chain_pdb),
        ignore_hydrogen=True,
        terminal_residues=[{"chain_id": "", "res_seq": 2, "n_terminal": True}],
        infer_terminals=False,
    )

    assert [res.name for res in blank_chain.residues] == ["NVAL"]


def test_pdb_ssbond_and_link_records_create_residue_links(tmp_path):
    _load_amber()
    pdb = """\
SSBOND   1 CYS A   1    CYS B   1
LINK         C   ACE C   1                 N   ALA C   2
ATOM      1  N   CYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  CYS A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   CYS A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   CYS A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  SG  CYS A   1       1.000   1.800   0.000  1.00  0.00           S
TER
ATOM      6  N   CYS B   1       0.000   5.000   0.000  1.00  0.00           N
ATOM      7  CA  CYS B   1       1.000   5.000   0.000  1.00  0.00           C
ATOM      8  C   CYS B   1       2.000   5.000   0.000  1.00  0.00           C
ATOM      9  O   CYS B   1       3.000   5.000   0.000  1.00  0.00           O
ATOM     10  SG  CYS B   1       1.000   3.200   0.000  1.00  0.00           S
TER
ATOM     11  CH3 ACE C   1       0.000  10.000   0.000  1.00  0.00           C
ATOM     12  C   ACE C   1       1.000  10.000   0.000  1.00  0.00           C
ATOM     13  O   ACE C   1       2.000  10.000   0.000  1.00  0.00           O
ATOM     14  N   ALA C   2       3.000  10.000   0.000  1.00  0.00           N
ATOM     15  CA  ALA C   2       4.000  10.000   0.000  1.00  0.00           C
ATOM     16  C   ALA C   2       5.000  10.000   0.000  1.00  0.00           C
ATOM     17  O   ALA C   2       6.000  10.000   0.000  1.00  0.00           O
TER
END
"""
    mol = Xponge.load_pdb(StringIO(pdb))

    Xponge.Save_SPONGE_Input(mol, prefix="linked", dirname=str(tmp_path))
    pairs = _bond_pairs(tmp_path / "linked_bond.txt")

    assert [res.name for res in mol.residues[:2]] == ["NCYX", "NCYX"]
    assert (4, 9) in pairs
    assert (11, 13) in pairs


def test_pdb_options_for_altloc_hydrogen_conect_and_cryst1(tmp_path):
    _load_amber()
    pdb = """\
CRYST1   10.000   11.000   12.000  90.00  91.00  92.00 P 1           1
ATOM      1  N  AALA A   1       0.000   0.000   0.000  0.50 10.00           N
ATOM      2  N  BALA A   1       9.000   9.000   9.000  0.50 20.00           N
ATOM      3  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      4  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      5  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      6  H   ALA A   1      -1.000   0.000   0.000  1.00  0.00           H
TER
HETATM    7  C1  LIG B   1       0.000   4.000   0.000  1.00  0.00           C
HETATM    8  C2  LIG B   2       1.000   4.000   0.000  1.00  0.00           C
CONECT    7    8
END
"""
    mol = Xponge.load_pdb(
        StringIO(pdb),
        position_need="B",
        ignore_hydrogen=True,
        ignore_conect=False,
        read_cryst1=True,
    )

    Xponge.Save_SPONGE_Input(mol, prefix="opts", dirname=str(tmp_path))

    assert mol.box_length == pytest.approx([10.0, 11.0, 12.0])
    assert mol.box_angle == pytest.approx([90.0, 91.0, 92.0])
    assert mol.residues[0].atoms[0].x == pytest.approx(9.0)
    assert all(atom.element != "H" for residue in mol.residues for atom in residue.atoms)
    assert (tmp_path / "opts_bond.txt").read_text().splitlines()[0] == "4"
    assert (4, 5) in _bond_pairs(tmp_path / "opts_bond.txt")


def test_pdb_writer_preserves_chain_boundaries_and_links(tmp_path):
    _load_amber()
    mol = Xponge.load_pdb(
        StringIO(
            """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       7.000   0.000   0.000  1.00  0.00           O
TER
ATOM      9  N   ALA B   1       0.000   5.000   0.000  1.00  0.00           N
ATOM     10  CA  ALA B   1       1.000   5.000   0.000  1.00  0.00           C
ATOM     11  C   ALA B   1       2.000   5.000   0.000  1.00  0.00           C
ATOM     12  O   ALA B   1       3.000   5.000   0.000  1.00  0.00           O
ATOM     13  N   GLY B   2       4.000   5.000   0.000  1.00  0.00           N
ATOM     14  CA  GLY B   2       5.000   5.000   0.000  1.00  0.00           C
ATOM     15  C   GLY B   2       6.000   5.000   0.000  1.00  0.00           C
ATOM     16  O   GLY B   2       7.000   5.000   0.000  1.00  0.00           O
TER
END
"""
        ),
        ignore_conect=False,
    )
    out = tmp_path / "roundtrip.pdb"

    Xponge.save_pdb(mol, str(out))
    text = out.read_text()
    reloaded = Xponge.load_pdb(str(out))

    assert text.startswith("REMARK   Generated By Xponge (Molecule)\n")
    assert text.count("\nTER\n") == 2
    assert " A   1" in text
    assert " B   1" in text
    assert [res.name for res in reloaded.residues] == ["NALA", "CGLY", "NALA", "CGLY"]


def test_pdb_writer_rebuilds_ssbond_link_and_conect_records(tmp_path):
    _load_amber()
    mol = Xponge.load_pdb(
        StringIO(
            """\
SSBOND   1 CYS A   1    CYS B   1
ATOM      1  N   CYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  CYS A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   CYS A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   CYS A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  SG  CYS A   1       1.000   1.800   0.000  1.00  0.00           S
ATOM      6  N   ALA A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      7  CA  ALA A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      8  C   ALA A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      9  O   ALA A   2       7.000   0.000   0.000  1.00  0.00           O
TER
ATOM     10  N   CYS B   1       0.000   5.000   0.000  1.00  0.00           N
ATOM     11  CA  CYS B   1       1.000   5.000   0.000  1.00  0.00           C
ATOM     12  C   CYS B   1       2.000   5.000   0.000  1.00  0.00           C
ATOM     13  O   CYS B   1       3.000   5.000   0.000  1.00  0.00           O
ATOM     14  SG  CYS B   1       1.000   3.200   0.000  1.00  0.00           S
ATOM     15  N   ALA B   2       4.000   5.000   0.000  1.00  0.00           N
ATOM     16  CA  ALA B   2       5.000   5.000   0.000  1.00  0.00           C
ATOM     17  C   ALA B   2       6.000   5.000   0.000  1.00  0.00           C
ATOM     18  O   ALA B   2       7.000   5.000   0.000  1.00  0.00           O
TER
CONECT    7   16
END
"""
        ),
        ignore_conect=False,
    )
    out = tmp_path / "links.pdb"

    Xponge.save_pdb(mol, str(out))
    text = out.read_text()
    reloaded = Xponge.load_pdb(str(out))

    assert "SSBOND   1 CYX A   1    CYX B   1" in text
    assert any(line.startswith("LINK") and "ALA A   2" in line and "ALA B   2" in line for line in text.splitlines())
    assert "CONECT" not in text
    Xponge.Save_SPONGE_Input(reloaded, prefix="roundtrip", dirname=str(tmp_path))
    assert (4, 13) in _bond_pairs(tmp_path / "roundtrip_bond.txt")
    assert (6, 15) in _bond_pairs(tmp_path / "roundtrip_bond.txt")


def test_pdb_writer_uses_conect_for_single_residue_chains(tmp_path):
    _load_amber()
    mol = Xponge.load_pdb(
        StringIO(
            """\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000  0.50 10.00           C
HETATM    2  C2  LIG B   1       1.000   0.000   0.000  0.75 20.00           C
CONECT    1    2
END
"""
        ),
        ignore_conect=False,
    )
    out = tmp_path / "single_link.pdb"

    Xponge.save_pdb(mol, str(out))
    text = out.read_text()
    reloaded = Xponge.load_pdb(str(out), ignore_conect=False)

    assert "LINK" not in text
    assert "SSBOND" not in text
    assert "CONECT    1    2" in text
    assert text.count("\nTER\n") == 2


def test_pdb_writer_hybrid36_resseq_and_link_records(tmp_path):
    mol2 = tmp_path / "h36.mol2"
    mol2.write_text(
        """@<TRIPOS>MOLECULE
H36
1 0 1 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A 0.0 0.0 0.0 C.3 1 H36 0.0
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
1 H36 1 **** 0 **** ****
"""
    )

    Xponge.load_mol2(str(mol2), as_template=True)
    Xponge.configure_residue_template_head("H36", "A")
    Xponge.configure_residue_template_tail("H36", "A")

    residue_type = Xponge.ResidueType.get_type("H36")
    mol = Xponge.Molecule("H36M")
    for _ in range(10001):
        mol.add_residue(Xponge.Residue(residue_type, directly_copy=True))
    for i in range(len(mol.residues) - 1):
        mol.add_residue_link(_residue_atom(mol.residues[i], "A"), _residue_atom(mol.residues[i + 1], "A"))
    mol.add_residue_link(_residue_atom(mol.residues[0], "A"), _residue_atom(mol.residues[-1], "A"))

    outfile = tmp_path / "h36m.pdb"
    Xponge.save_pdb(mol, str(outfile))
    lines = outfile.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    link_lines = [line for line in lines if line.startswith("LINK")]

    assert atom_lines[9998][22:26] == "9999"
    assert atom_lines[9999][22:26] == "A000"
    assert atom_lines[10000][22:26] == "A001"
    assert len(link_lines) == 1
    assert link_lines[0][22:26] == "   1"
    assert link_lines[0][52:56] == "A001"

    roundtrip = Xponge.load_pdb(str(outfile), ignore_conect=False)
    assert len(roundtrip.residues) == 10001


def test_assign_save_as_pdb_uses_hybrid36_serials_and_conect(tmp_path):
    assign = Xponge.Assign("BEN")
    for _ in range(100001):
        assign.add_atom("C", 0.0, 0.0, 0.0, name="C")
    assign.add_bond(99999, 100000, 1)

    outfile = tmp_path / "assign_h36.pdb"
    assign.save_as_pdb(str(outfile))
    lines = outfile.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    conect_lines = [line for line in lines if line.startswith("CONECT")]

    assert atom_lines[99999][6:11] == "A0000"
    assert atom_lines[100000][6:11] == "A0001"
    assert any(line[6:11] == "A0000" and line[11:16] == "A0001" for line in conect_lines)


def test_pdb_writer_outputs_seqres_only_for_multi_residue_chains(tmp_path):
    _load_amber()
    mol = Xponge.load_pdb(
        StringIO(
            """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       7.000   0.000   0.000  1.00  0.00           O
TER
HETATM    9  C1  LIG Z   1       0.000   4.000   0.000  1.00  0.00           C
END
"""
        )
    )
    out = tmp_path / "seqres.pdb"

    Xponge.save_pdb(mol, str(out))
    seqres = [line for line in out.read_text().splitlines() if line.startswith("SEQRES")]

    assert seqres == ["SEQRES   1 A    2  ALA GLY"]


def test_pdb_hybrid36_read_large_indices():
    _load_amber()
    decoded = Xponge.load_pdb(
        StringIO(
            "ATOM  A0000  N   ALA AA000       0.000   0.000   0.000  1.00  0.00           N\n"
            "END\n"
        ),
        unterminal_residues=["A:10000"],
    )
    assert decoded.residues[0].pdb_resseq == 10000
    assert decoded.residues[0].atoms[0].serial == 100000


def test_molecule_plus_and_pipe_follow_xponge_link_semantics(tmp_path):
    _load_amber()
    ala = Xponge.get_template_molecule("ALA")
    gly = Xponge.get_template_molecule("GLY")

    linked = ala + gly
    unlinked = ala | gly
    repeated = ala * 3

    Xponge.Save_SPONGE_Input(linked, prefix="linked", dirname=str(tmp_path))
    Xponge.Save_SPONGE_Input(unlinked, prefix="unlinked", dirname=str(tmp_path))
    Xponge.Save_SPONGE_Input(repeated, prefix="repeated", dirname=str(tmp_path))

    linked_pairs = _bond_pairs(tmp_path / "linked_bond.txt")
    unlinked_pairs = _bond_pairs(tmp_path / "unlinked_bond.txt")

    assert linked.residue_count == 2
    assert unlinked.residue_count == 2
    assert repeated.residue_count == 3
    assert any(pair[0] < ala.atom_count <= pair[1] for pair in linked_pairs)
    assert not any(pair[0] < ala.atom_count <= pair[1] for pair in unlinked_pairs)
