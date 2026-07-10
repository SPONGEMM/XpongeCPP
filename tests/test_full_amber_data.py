import importlib
from importlib import resources
from io import StringIO
from pathlib import Path

import pytest

import XpongeCPP as Xponge


def test_full_amber_forcefield_data_is_packaged():
    amber_root = resources.files("XpongeCPP").joinpath("data", "amber")

    required = [
        "ff14SB.mol2",
        "ff14SB.frcmod",
        "parm10.dat",
        "ff19SB.mol2",
        "ff19SB.frcmod",
        "gaff.dat",
        "gaff2.dat",
        "tip3p.mol2",
        "atomic_ions.mol2",
        "ions1lm_126_tip3p.frcmod",
        "ionsjc_tip3p.frcmod",
        "ions234lm_126_tip3p.frcmod",
        "opc.mol2",
        "spce.mol2",
        "glycam_06j/GLYCAM_06j.dat",
    ]

    for relative in required:
        assert amber_root.joinpath(relative).is_file(), relative


def test_reference_xponge_forcefield_tree_is_available_in_source_tree():
    reference_root = Path(__file__).resolve().parents[1] / "third_party" / "xponge_reference_forcefield"

    required = [
        "amber/ff14SB.mol2",
        "amber/gaff2.dat",
        "charmm/charmm36/forcefield.itp",
        "opls/oplsaam/forcefield.itp",
        "martini/martini300/martini_v3.0.0.itp",
    ]

    for relative in required:
        assert reference_root.joinpath(relative).is_file(), relative


def test_ff14sb_and_tip3p_imports_register_templates_from_real_amber_mol2():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    assert Xponge.template_atom_count("MET") >= 17
    assert Xponge.template_atom_count("SER") >= 11
    assert Xponge.template_atom_count("WAT") == 3
    assert Xponge.template_atom_count("NA") == 1
    assert Xponge.template_atom_count("CL") == 1

    water = Xponge.get_template_molecule("WAT")
    assert water.atom_count == 3
    assert water.residues[0].name == "WAT"
    assert water.residues[0].name2atom("O").type == "OW"


def test_ff19sb_import_registers_real_templates_and_cmap_parameters(tmp_path):
    import XpongeCPP.forcefield.amber.ff19sb  # noqa: F401

    assert Xponge.template_atom_count("ALA") >= 10
    assert Xponge.template_atom_count("GLY") >= 7

    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
FF19_CMAP_TEST
5 4 3
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C  0.0 0.0 0.0 C 1 PRE 0.0
2 N  1.3 0.0 0.0 N 2 ALA 0.0
3 CA 2.6 0.0 0.0 XC 2 ALA 0.0
4 C  3.9 0.0 0.0 C 2 ALA 0.0
5 N  5.2 0.0 0.0 N 3 NXT 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
4 4 5 1
"""
        )
    )

    out = Xponge.Save_SPONGE_Input(mol, prefix="ff19", dirname=str(tmp_path))

    assert "cmap" in out
    assert (tmp_path / "ff19_cmap.txt").read_text().splitlines()[0] == "1 1"


def test_gaff_and_gaff2_imports_register_packaged_parameters(tmp_path):
    import XpongeCPP.forcefield.amber.gaff  # noqa: F401
    import XpongeCPP.forcefield.amber.gaff2  # noqa: F401

    # A missing parameter would raise during export; this exercises the Python import entry points.
    mol = Xponge.load_mol2(
        StringIO(
        """@<TRIPOS>MOLECULE
ETH
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 c3 1 ETH 0.0
2 C2 1.5 0.0 0.0 c3 1 ETH 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    assert mol.atom_count == 2
    assert mol.validate()
    Xponge.Save_SPONGE_Input(mol, prefix="eth", dirname=str(tmp_path))
    assert (tmp_path / "eth_bond.txt").read_text().splitlines()[0] == "1"


def test_amber_multi_site_water_imports_register_complete_templates(tmp_path):
    for module_name, prefix, expected_k in [
        ("XpongeCPP.forcefield.amber.tip4p", "tip4p", "0.127970"),
        ("XpongeCPP.forcefield.amber.tip4pew", "tip4pew", "0.106641"),
        ("XpongeCPP.forcefield.amber.opc", "opc", "0.147721"),
    ]:
        module = importlib.import_module(module_name)
        importlib.reload(module)

        water = Xponge.get_template_molecule("WAT")
        assert water.atom_count == 4
        assert water.residues[0].name2atom("EPW").type == "EP"

        out = Xponge.Save_SPONGE_Input(water, prefix=prefix, dirname=str(tmp_path))
        assert "virtual_atom" in out
        assert (tmp_path / f"{prefix}_virtual_atom.txt").read_text().splitlines() == [
            f"2 3 0 1 2 {expected_k} {expected_k}",
        ]


def test_spce_import_registers_three_site_water_without_virtual_atom(tmp_path):
    import XpongeCPP.forcefield.amber.spce  # noqa: F401

    water = Xponge.get_template_molecule("WAT")
    assert water.atom_count == 3

    out = Xponge.Save_SPONGE_Input(water, prefix="spce", dirname=str(tmp_path))
    assert "virtual_atom" not in out
    assert (tmp_path / "spce_bond.txt").read_text().splitlines()[0] == "3"


def test_water_model_ion_replacement_workflows_validate_across_supported_amber_models():
    solute = Xponge.load_pdb(
        StringIO(
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        )
    )

    for module_name, expected_atoms_per_water in [
        ("XpongeCPP.forcefield.amber.tip3p", 3),
        ("XpongeCPP.forcefield.amber.spce", 3),
        ("XpongeCPP.forcefield.amber.tip4p", 4),
        ("XpongeCPP.forcefield.amber.tip4pew", 4),
        ("XpongeCPP.forcefield.amber.opc", 4),
    ]:
        module = importlib.import_module(module_name)
        importlib.reload(module)

        mol = solute.copy()
        water = Xponge.get_template_molecule("WAT")
        assert water.atom_count == expected_atoms_per_water

        Xponge.Add_Solvent_Box(mol, water, 4.0, tolerance=2.5, n_solvent=4, seed=20260509)
        Xponge.Add_Ions(mol, {"NA": 1}, seed=20260509)

        counts = mol.residue_counts()
        assert counts["NA"] == 1
        assert counts["WAT"] == 3
        assert mol.validate()


def test_amber_frcmod_cmap_is_generated_automatically(tmp_path):
    frcmod = tmp_path / "mini_cmap.frcmod"
    frcmod.write_text(
        """mini cmap
CMAP
%FLAG CMAP_COUNT   1
%FLAG CMAP_RESLIST  1
ALA
%FLAG CMAP_RESOLUTION 2
%FLAG CMAP_PARAMETER
 1.0 2.0
 3.0 4.0
"""
    )
    Xponge.register_amber_frcmod_file(str(frcmod))
    for atom_type, rmin, epsilon in [("C", 1.7, 0.1), ("N", 1.5, 0.1), ("XC", 1.8, 0.1)]:
        Xponge.register_amber_lj_parameter(atom_type, atom_type, epsilon, rmin)

    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
CMAP_TEST
5 4 3
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C  0.0 0.0 0.0 C 1 PRE 0.0
2 N  1.3 0.0 0.0 N 2 ALA 0.0
3 CA 2.6 0.0 0.0 XC 2 ALA 0.0
4 C  3.9 0.0 0.0 C 2 ALA 0.0
5 N  5.2 0.0 0.0 N 3 NXT 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
4 4 5 1
"""
        )
    )

    out = Xponge.Save_SPONGE_Input(mol, prefix="cmap", dirname=str(tmp_path))

    assert "cmap" in out
    assert (tmp_path / "cmap_cmap.txt").read_text().splitlines() == [
        "1 1",
        "2 ",
        "1.000000 2.000000 ",
        "3.000000 4.000000 ",
        "",
        "0 1 2 3 4 0",
    ]


def test_amber_frcmod_dihedral_continuation_and_nb14(tmp_path):
    frcmod = tmp_path / "continuation.frcmod"
    frcmod.write_text(
        """continuation test
DIHE
qA-qB-qC-qD    1  0.6  0.0  -3.0  SCEE=1.0 SCNB=1.0
              1  0.2  180.0  1.0  SCEE=1.0 SCNB=1.0

NONB
qA  1.5  0.1
qB  1.5  0.1
qC  1.5  0.1
qD  1.5  0.1
"""
    )
    Xponge.register_amber_frcmod_file(str(frcmod))
    for left, right in [("qA", "qB"), ("qB", "qC"), ("qC", "qD")]:
        Xponge.register_amber_bond_parameter(left, right, 100.0, 1.5)
    Xponge.register_amber_angle_parameter(["qA", "qB", "qC"], 50.0, 2.0)
    Xponge.register_amber_angle_parameter(["qB", "qC", "qD"], 50.0, 2.0)

    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
CONTINUATION
4 3 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A 0.0 0.0 0.0 qA 1 MOL 0.1
2 B 1.5 0.0 0.0 qB 1 MOL -0.1
3 C 3.0 0.0 0.0 qC 1 MOL 0.1
4 D 4.5 0.0 0.0 qD 1 MOL -0.1
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
"""
        )
    )
    Xponge.Save_SPONGE_Input(mol, prefix="continuation", dirname=str(tmp_path))

    dihedrals = (tmp_path / "continuation_dihedral.txt").read_text().splitlines()
    nb14 = (tmp_path / "continuation_nb14.txt").read_text().splitlines()
    assert dihedrals[0] == "2"
    assert {int(line.split()[4]) for line in dihedrals[1:]} == {1, 3}
    assert nb14[1].split()[2:] == ["1.000000", "1.000000"]


def test_amber_frcmod_rejects_orphan_dihedral_continuation(tmp_path):
    frcmod = tmp_path / "orphan.frcmod"
    frcmod.write_text(
        """orphan continuation test
DIHE
              1  0.2  180.0  1.0
"""
    )

    with pytest.raises(RuntimeError, match="no preceding atom types"):
        Xponge.register_amber_frcmod_file(str(frcmod))


def test_full_amber_ecosystem_import_modules_register_packaged_templates_and_parameters():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    for module_name in [
        "XpongeCPP.forcefield.amber.bsc1",
        "XpongeCPP.forcefield.amber.ol3",
        "XpongeCPP.forcefield.amber.ol15",
        "XpongeCPP.forcefield.amber.lipid14",
        "XpongeCPP.forcefield.amber.lipid17",
        "XpongeCPP.forcefield.amber.rsff2c",
        "XpongeCPP.forcefield.amber.glycam_06j.d_furanose",
        "XpongeCPP.forcefield.amber.glycam_06j.d_pyranose",
        "XpongeCPP.forcefield.amber.glycam_06j.l_furanose",
        "XpongeCPP.forcefield.amber.glycam_06j.l_pyranose",
        "XpongeCPP.forcefield.amber.glycam_06j.glycoprotein",
    ]:
        importlib.import_module(module_name)

    assert Xponge.template_atom_count("DA") > 20
    assert Xponge.template_atom_count("A") > 20
    assert Xponge.template_atom_count("ACE") == 6
    assert Xponge.template_atom_count("NME") == 6


def test_amber_nucleic_lipid_and_glycam_modules_support_representative_export_workflows(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.bsc1  # noqa: F401
    import XpongeCPP.forcefield.amber.ol3  # noqa: F401
    import XpongeCPP.forcefield.amber.lipid17  # noqa: F401
    import XpongeCPP.forcefield.amber.glycam_06j.d_pyranose  # noqa: F401

    cases = [
        ("dna", Xponge.get_template_molecule("DA5") + Xponge.get_template_molecule("DT3"), ["DA5", "DT3"]),
        ("rna", Xponge.get_template_molecule("A5") + Xponge.get_template_molecule("U3"), ["A5", "U3"]),
        ("lipid", Xponge.get_template_molecule("PA"), ["PA"]),
        ("glycam", Xponge.get_template_molecule("0GA") + Xponge.get_template_molecule("4YB"), ["0GA", "4YB"]),
    ]

    for prefix, mol, expected_resnames in cases:
        assert [res.name for res in mol.residues] == expected_resnames
        assert mol.validate()

        out = Xponge.Save_SPONGE_Input(mol, prefix=prefix, dirname=str(tmp_path))

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
        assert int((tmp_path / f"{prefix}_bond.txt").read_text().splitlines()[0]) > 0
        assert (tmp_path / f"{prefix}_resname.txt").read_text().splitlines() == [str(len(expected_resnames)), *expected_resnames]


def test_amber_nucleic_lipid_and_glycam_modules_support_broader_assembled_export_workflows(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.bsc1  # noqa: F401
    import XpongeCPP.forcefield.amber.ol3  # noqa: F401
    import XpongeCPP.forcefield.amber.lipid17  # noqa: F401
    import XpongeCPP.forcefield.amber.glycam_06j.d_pyranose  # noqa: F401

    cases = [
        ("dna3", Xponge.get_template_molecule("DA5") + Xponge.get_template_molecule("DT") + Xponge.get_template_molecule("DG3"), ["DA5", "DT", "DG3"]),
        ("rna3", Xponge.get_template_molecule("A5") + Xponge.get_template_molecule("U") + Xponge.get_template_molecule("G3"), ["A5", "U", "G3"]),
        ("lipid2", Xponge.get_template_molecule("PC") + Xponge.get_template_molecule("OL"), ["PC", "OL"]),
        ("glycam3", Xponge.get_template_molecule("0GA") + Xponge.get_template_molecule("4YB") + Xponge.get_template_molecule("4YB"), ["0GA", "4YB", "4YB"]),
    ]

    for prefix, mol, expected_resnames in cases:
        assert [res.name for res in mol.residues] == expected_resnames
        assert mol.atom_count > len(expected_resnames)
        assert mol.validate()

        out = Xponge.Save_SPONGE_Input(mol, prefix=prefix, dirname=str(tmp_path))

        assert {"bond", "angle", "dihedral", "exclude", "nb14", "residue", "resname"}.issubset(out)
        assert int((tmp_path / f"{prefix}_bond.txt").read_text().splitlines()[0]) > 0
        assert int((tmp_path / f"{prefix}_angle.txt").read_text().splitlines()[0]) > 0
        assert int((tmp_path / f"{prefix}_dihedral.txt").read_text().splitlines()[0]) > 0
        assert (tmp_path / f"{prefix}_resname.txt").read_text().splitlines() == [str(len(expected_resnames)), *expected_resnames]


def test_pdb_terminal_mapping_is_chain_local_and_keeps_ace_nme_caps():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    pdb = """\
ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  MET A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   MET A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   MET A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   SER A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  SER A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   SER A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   SER A   2       7.000   0.000   0.000  1.00  0.00           O
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
ATOM     17  CH3 ACE C   1       0.000  10.000   0.000  1.00  0.00           C
ATOM     18  C   ACE C   1       1.000  10.000   0.000  1.00  0.00           C
ATOM     19  O   ACE C   1       2.000  10.000   0.000  1.00  0.00           O
ATOM     20  N   ALA C   2       3.000  10.000   0.000  1.00  0.00           N
ATOM     21  CA  ALA C   2       4.000  10.000   0.000  1.00  0.00           C
ATOM     22  C   ALA C   2       5.000  10.000   0.000  1.00  0.00           C
ATOM     23  O   ALA C   2       6.000  10.000   0.000  1.00  0.00           O
ATOM     24  N   NME C   3       7.000  10.000   0.000  1.00  0.00           N
ATOM     25  CH3 NME C   3       8.000  10.000   0.000  1.00  0.00           C
TER
END
"""

    mol = Xponge.load_pdb(StringIO(pdb))

    assert [res.name for res in mol.residues] == [
        "NMET",
        "CSER",
        "NALA",
        "CGLY",
        "ACE",
        "ALA",
        "NME",
    ]
