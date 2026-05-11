from io import StringIO

import XpongeCPP as Xponge


def test_gromacs_topology_parser_generates_special_forces(tmp_path):
    include = tmp_path / "params.itp"
    include.write_text(
        """[ atomtypes ]
; name bond_type atomic_number mass charge ptype sigma epsilon
Q Q 6 12.011 0.0 A 3.50000E-01 2.76144E-01
H H 1 1.008 0.0 A 2.50000E-01 1.25520E-01
"""
    )
    top = tmp_path / "system.top"
    top.write_text(
        f"""#include "{include.name}"
[ moleculetype ]
MOL 3
[ atoms ]
1 Q 1 MOL A1 1 0.0 12.011
2 Q 1 MOL A2 1 0.0 12.011
3 Q 1 MOL A3 1 0.0 12.011
4 H 1 MOL A4 1 0.0 1.008
5 H 1 MOL VS 1 0.0 0.0
[ bonds ]
1 2 1
2 3 1
3 4 1
[ pairs ]
1 4 1 1.25 2.5 0.75
[ angles ]
1 2 3 5 109.5 50.0 2.0 10.0
[ dihedrals ]
1 2 3 4 3 0.1 0.2 0.3 0.4 0.5 0.6
[ virtual_sites3 ]
5 1 2 3 1 0.25 0.75
"""
    )

    mol = Xponge.load_gromacs_topology_file(str(top))
    out = Xponge.Save_SPONGE_Input(mol, prefix="gmx", dirname=str(tmp_path))

    assert {"urey_bradley", "Ryckaert_Bellemans", "nb14_extra", "virtual_atom"}.issubset(out)
    assert (tmp_path / "gmx_urey_bradley.txt").read_text().splitlines()[0] == "1"
    assert (tmp_path / "gmx_Ryckaert_Bellemans.txt").read_text().splitlines()[0] == "1"
    assert (tmp_path / "gmx_nb14_extra.txt").read_text().splitlines()[0] == "1"
    assert (tmp_path / "gmx_virtual_atom.txt").read_text().splitlines() == [
        "2 4 0 1 2 0.250000 0.750000",
    ]


def test_opls_itp_parser_generates_rb_dihedral(tmp_path):
    itp = tmp_path / "opls.itp"
    itp.write_text(
        """[ atomtypes ]
opls_1 C 6 12.011 0.0 A 0.35 0.276144
[ moleculetype ]
MOL 3
[ atoms ]
1 opls_1 1 MOL A1 1 0.0 12.011
2 opls_1 1 MOL A2 1 0.0 12.011
3 opls_1 1 MOL A3 1 0.0 12.011
4 opls_1 1 MOL A4 1 0.0 12.011
[ bonds ]
1 2 1
2 3 1
3 4 1
[ dihedrals ]
1 2 3 4 3 0.1 0.2 0.3 0.4 0.5 0.6
"""
    )

    mol = Xponge.load_opls_itp_file(str(itp))
    out = Xponge.Save_SPONGE_Input(mol, prefix="opls", dirname=str(tmp_path))

    assert "Ryckaert_Bellemans" in out
    assert (tmp_path / "opls_Ryckaert_Bellemans.txt").read_text().splitlines() == [
        "1",
        "0 1 2 3 0.100000 0.200000 0.300000 0.400000 0.500000 0.600000",
    ]


def test_charmm_prm_and_rtf_parser_generates_template_and_nbfix(tmp_path):
    prm = tmp_path / "mini.prm"
    prm.write_text(
        """MASS -1 CT 12.011
MASS -1 HT 1.008
NONBONDED
CT 0.0 -0.100 2.000
HT 0.0 -0.050 1.000
ANGLE
HT CT HT 50.0 109.5 10.0 2.0
NBFIX
CT HT 1.25 2.5 0.75
"""
    )
    rtf = tmp_path / "mini.rtf"
    rtf.write_text(
        """RESI ETH 0.0
ATOM H1 HT 0.1
ATOM C1 CT -0.2
ATOM H2 HT 0.1
BOND H1 C1 C1 H2
END
"""
    )

    Xponge.load_charmm_parameter_file(str(prm))
    mol = Xponge.load_charmm_topology_file(str(rtf))
    out = Xponge.Save_SPONGE_Input(mol, prefix="charmm", dirname=str(tmp_path))

    assert mol.residue_count == 1
    assert mol.atom_count == 3
    assert {"urey_bradley", "nb14_extra"}.issubset(out)
    assert (tmp_path / "charmm_urey_bradley.txt").read_text().splitlines()[0] == "1"
    assert (tmp_path / "charmm_nb14_extra.txt").read_text().splitlines()[0] == "2"


def test_sw_and_edip_parameter_parsers_bind_atom_types(tmp_path):
    mol = Xponge.load_mol2(
        StringIO(
        """@<TRIPOS>MOLECULE
PAIRWISE
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A 0.0 0.0 0.0 S 1 SOL 0.0
2 B 1.0 0.0 0.0 S 1 SOL 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    sw = tmp_path / "mw.sw"
    sw.write_text(
        """S-S 1.1 2.2 3.3 4.0 5.0 6.6 7.7 8.8 0.0 0.0
S-S-S 0.0 0.0 3.3 0.0 0.0 0.0 0.0 0.0 9.9 10.1
"""
    )
    edip = tmp_path / "si.edip"
    edip.write_text(
        """S-S 1.1 2.2 3.3 4.4 5.5 6.6 7.7 8.8 0.0 0.0 0.0 12.12 0.0 0.0 0.0 0.0 0.0
S-S-S 0.0 0.0 0.0 0.0 0.0 0.0 9.9 10.1 11.11 12.12 13.13 0.0 14.14 15.15 16.16 17.17 18.18
"""
    )

    Xponge.load_sw_parameter_file(str(sw), mol)
    Xponge.load_edip_parameter_file(str(edip), mol)
    out = Xponge.Save_SPONGE_Input(mol, prefix="pair", dirname=str(tmp_path))

    assert {"SW", "EDIP"}.issubset(out)
    assert (tmp_path / "pair_SW.txt").read_text().splitlines()[0] == "2 1"
    assert (tmp_path / "pair_EDIP.txt").read_text().splitlines()[0] == "2 1"


def test_softcore_export_uses_b_types_and_subsystems(tmp_path):
    Xponge.register_amber_lj_parameter("A", "A", 0.2, 1.0)
    Xponge.register_amber_lj_parameter("B", "B", 0.3, 1.5)
    mol = Xponge.load_mol2(
        StringIO(
        """@<TRIPOS>MOLECULE
FEP
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 A 1 FEP 0.0
2 A2 1.0 0.0 0.0 A 1 FEP 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )
    for atom in mol.residues[0].atoms:
        atom.lj_type_b = "B"
        atom.subsys = 1
    mol.enable_lj_soft_core()

    out = Xponge.Save_SPONGE_Input(mol, prefix="fep", dirname=str(tmp_path))

    assert {"LJ_soft_core", "subsys_division"}.issubset(out)
    assert "LJ" not in out
    assert (tmp_path / "fep_subsys_division.txt").read_text().splitlines() == ["2", "1", "1"]


def test_non_amber_forcefield_import_modules_are_available():
    import XpongeCPP.forcefield.charmm.charmm36  # noqa: F401
    import XpongeCPP.forcefield.opls.oplsaam  # noqa: F401
    import XpongeCPP.forcefield.sw.mw as mw
    import XpongeCPP.forcefield.edip.si as si

    assert callable(mw.load_parameters)
    assert callable(si.load_parameters)
