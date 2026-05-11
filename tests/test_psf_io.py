from io import StringIO

import XpongeCPP as Xponge


def test_load_molpsf_splits_residues_and_molecules_by_connectivity(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    psf = StringIO(
        """\
PSF

       1 !NTITLE
 REMARKS fixture
       4 !NATOM
       1 SEG1     1 WAT  O1   OW     -0.834000       15.9990           0
       2 SEG1     1 WAT  H1   HW      0.417000        1.0080           0
       3 SEG1     1 WAT  O2   OW     -0.834000       15.9990           0
       4 SEG1     1 WAT  H2   HW      0.417000        1.0080           0

       2 !NBOND: bonds
       1       2       3       4
"""
    )

    mol, mols = Xponge.load_molpsf(psf)

    assert mol.atom_count == 4
    assert mol.residue_count == 2
    assert sorted(mols) == ["psf_1", "psf_2"]
    assert mols["psf_1"].atom_count == 2
    assert mols["psf_2"].atom_count == 2
    Xponge.Save_SPONGE_Input(mol, prefix="psf", dirname=str(tmp_path))
    assert (tmp_path / "psf_bond.txt").read_text().splitlines()[0] == "2"


def test_load_molpsf_can_keep_single_system_without_connectivity_split():
    psf = StringIO(
        """\
PSF
       0 !NTITLE
       2 !NATOM
       1 SYS      1 LIG  C1   CT      0.000000       12.0100           0
       2 SYS      2 LIG  C2   CT      0.000000       12.0100           0
       1 !NBOND: bonds
       1       2
"""
    )

    mol, mols = Xponge.load_molpsf(psf, split_by=None)

    assert mol.atom_count == 2
    assert mol.residue_count == 2
    assert list(mols) == ["psf"]


def test_load_molpsf_accepts_header_variants_and_skips_unused_sections(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    psf = StringIO(
        """\
PSF EXT XPLOR CHARMM
       0 !NTITLE
       3 !NATOM
       1 SYS      1 WAT  O    OW     -0.834000       15.9990           0
       2 SYS      1 WAT  H1   HW      0.417000        1.0080           0
       3 SYS      1 WAT  H2   HW      0.417000        1.0080           0
       2 !NBOND: bonds
       1       2       2       3
       1 !NTHETA: angles
       1       2       3
       1 !NPHI: dihedrals
       1       2       3       1
"""
    )

    mol, mols = Xponge.load_molpsf(psf)

    assert mol.atom_count == 3
    assert mol.residue_count == 1
    assert sorted(mols) == ["psf_1"]
    assert [atom.charge for atom in mol.residues[0].atoms] == [-0.834, 0.417, 0.417]
    Xponge.Save_SPONGE_Input(mol, prefix="psf_skip", dirname=str(tmp_path))
    assert (tmp_path / "psf_skip_bond.txt").read_text().splitlines()[0] == "2"
    assert (tmp_path / "psf_skip_angle.txt").read_text().splitlines()[0] == "0"


def test_load_molpsf_splits_interleaved_same_residue_components():
    psf = StringIO(
        """\
PSF
       0 !NTITLE
       4 !NATOM
       1 SYS      1 WAT  O1   OW     -0.834000       15.9990           0
       2 SYS      1 WAT  O2   OW     -0.834000       15.9990           0
       3 SYS      1 WAT  H1   HW      0.417000        1.0080           0
       4 SYS      1 WAT  H2   HW      0.417000        1.0080           0
       2 !NBOND: bonds
       1       3       2       4
"""
    )

    mol, mols = Xponge.load_molpsf(psf)

    assert mol.residue_count == 2
    assert [[atom.name for atom in residue.atoms] for residue in mol.residues] == [["O1", "H1"], ["O2", "H2"]]
    assert sorted(mols) == ["psf_1", "psf_2"]


def test_load_molpsf_distinguishes_reused_residue_name_with_charge_type_conflicts():
    psf = StringIO(
        """\
PSF
       0 !NTITLE
       2 !NATOM
       1 SYS      1 LIG  C1   CT      0.000000       12.0100           0
       2 SYS      2 LIG  C1   C2      0.500000       12.0100           0
       0 !NBOND: bonds
"""
    )

    mol, _ = Xponge.load_molpsf(psf, split_by=None)

    assert [res.type_name for res in mol.residues] == ["LIG", "LIG_1"]
    assert [res.atoms[0].type for res in mol.residues] == ["CT", "C2"]
    assert [res.atoms[0].charge for res in mol.residues] == [0.0, 0.5]
