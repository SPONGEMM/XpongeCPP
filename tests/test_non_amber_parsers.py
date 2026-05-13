import json
import subprocess
import sys
import textwrap
from io import StringIO

import XpongeCPP as Xponge
import pytest
from conftest import original_xponge_repo


XPONGE_REPO = original_xponge_repo()


def _run_original_xponge(script):
    if not XPONGE_REPO.exists():
        pytest.skip("local Xponge reference repository is not available")
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=XPONGE_REPO,
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        pytest.skip(f"local Xponge reference failed: {result.stderr[-500:]}")
    return result.stdout


def _pdb_atom_line(serial, name, resname, chain, resseq, x, y, z, element):
    return (
        f"{'ATOM':<6}{serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}"
        f"    {x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}          {element:>2}\n"
    )


def test_gromacs_topology_iterator_handles_macros_include_and_continuation(tmp_path):
    include = tmp_path / "included.itp"
    include.write_text(
        """#define VALUE 0.25
[ atomtypes ]
Q Q 6 12.011 0.0 A 0.35 \\
VALUE
#ifdef EXTRA
H H 1 1.008 0.0 A 0.25 0.125
#else
X X 0 0.0 0.0 A 0.10 0.010
#endif
"""
    )
    top = tmp_path / "topol.top"
    top.write_text(f"""#include "{include.name}"\n""")

    iterator = Xponge.GromacsTopologyIterator(str(top), macros={"EXTRA": ""})
    lines = list(iterator)

    assert iterator.defined_macros["VALUE"] == "0.25"
    assert lines == [
        "Q Q 6 12.011 0.0 A 0.35  0.25",
        "H H 1 1.008 0.0 A 0.25 0.125",
    ]


def test_load_ffitp_returns_xponge_style_forcefield_buffers(tmp_path):
    ffitp = tmp_path / "ffbonded.itp"
    ffitp.write_text(
        """[ defaults ]
1 2 yes 0.5 0.833333
[ atomtypes ]
Q Q 6 12.011 0.0 A 0.35 0.25
[ pairtypes ]
Q H 1 0.12 0.34
[ bondtypes ]
Q H 1 0.109 300000
[ angletypes ]
H Q H 5 109.5 50.0 0.2 10.0
[ dihedraltypes ]
Q Q Q Q 3 0.1 0.2 0.3 0.4 0.5 0.6
[ cmaptypes ]
Q Q Q Q Q 1 24 2 0.0 1.0 2.0 3.0
"""
    )

    output = Xponge.load_ffitp(str(ffitp))

    assert "Q 12.011 0.0 Q" in output["atomtypes"]
    assert "Q-Q 0.35 0.25" in output["LJ"]
    assert "Q-H 0.12 0.34 0.833333" in output["nb14_extra"]
    assert "Q-H 0.109 150000.0" in output["bonds"]
    assert "H-Q-H 109.5 25.0 0.2 5.0" in output["Urey-Bradley"]
    assert "Q-Q-Q-Q 0.1 0.2 0.3 0.4 0.5 0.6" in output["RB_dihedrals"]
    assert output["cmaps"]["Q-Q-Q-Q-Q"]["resolution"] == 2


def test_load_ffitp_matches_original_xponge_reference(tmp_path):
    ffitp = tmp_path / "ffbonded.itp"
    ffitp.write_text(
        """[ defaults ]
1 2 yes 0.5 0.833333
[ atomtypes ]
Q Q 6 12.011 0.0 A 0.35 0.25
H H 1 1.008 0.0 A 0.25 0.125
[ pairtypes ]
Q H 1 0.12 0.34
[ bondtypes ]
Q H 1 0.109 300000
[ angletypes ]
H Q H 5 109.5 50.0 0.2 10.0
[ dihedraltypes ]
Q Q Q Q 3 0.1 0.2 0.3 0.4 0.5 0.6
[ cmaptypes ]
Q Q Q Q Q 1 24 2 0.0 1.0 2.0 3.0
[ nonbond_params ]
Q H 1 0.20 0.10
"""
    )

    current = Xponge.load_ffitp(str(ffitp))
    script = textwrap.dedent(
        f"""
        import json
        import sys
        sys.path.insert(0, {str(XPONGE_REPO)!r})
        import Xponge

        output = Xponge.load_ffitp({str(ffitp)!r})
        payload = {{}}
        for key, value in output.items():
            if key == "cmaps":
                payload[key] = {{
                    subkey: {{
                        "resolution": int(subvalue["resolution"]),
                        "parameters": [float(item) for item in subvalue["parameters"]],
                    }}
                    for subkey, subvalue in value.items()
                }}
            elif key == "bond_type_names":
                payload[key] = dict(value)
            else:
                payload[key] = value
        print(json.dumps(payload, sort_keys=True))
        """
    )
    reference = json.loads(_run_original_xponge(script))

    current_payload = {}
    for key, value in current.items():
        if key == "cmaps":
            current_payload[key] = {
                subkey: {
                    "resolution": int(subvalue["resolution"]),
                    "parameters": [float(item) for item in subvalue["parameters"]],
                }
                for subkey, subvalue in value.items()
            }
        elif key == "bond_type_names":
            current_payload[key] = dict(value)
        else:
            current_payload[key] = value
    assert current_payload == reference


def test_load_molitp_builds_system_with_head_tail_residues_and_molecule_counts(tmp_path):
    top = tmp_path / "system.top"
    top.write_text(
        """[ moleculetype ]
PEP 3
[ atoms ]
1 N 1 ALA N 1 -0.3 14.0
2 C 1 ALA CA 1 0.1 12.0
3 C 2 GLY C 2 0.2 12.0
4 O 2 GLY O 2 -0.2 16.0
[ bonds ]
1 2 1
2 3 1
3 4 1
[ system ]
case
[ molecules ]
PEP 2
"""
    )

    system, mols = Xponge.load_molitp(str(top), water_replace=False)

    assert sorted(mols) == ["PEP"]
    assert mols["PEP"].residue_count == 2
    residue_names = [res.name for res in mols["PEP"].residues]
    assert residue_names[0].startswith("NALA")
    assert residue_names[1].startswith("CGLY")
    assert system.name == "case"
    assert system.atom_count == 8
    counts = system.residue_counts()
    assert counts[residue_names[0]] == 2
    assert counts[residue_names[1]] == 2


def test_load_molitp_matches_original_xponge_reference_for_multi_molecule_system(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    top = tmp_path / "system.top"
    top.write_text(
        """[ moleculetype ]
PEP 3
[ atoms ]
1 N 1 ALA N 1 -0.3 14.0
2 C 1 ALA CA 1 0.1 12.0
3 C 2 GLY C 2 0.2 12.0
4 O 2 GLY O 2 -0.2 16.0
[ bonds ]
1 2 1
2 3 1
3 4 1
[ moleculetype ]
SOL 3
[ atoms ]
1 OW 1 SOL OW 1 -0.834 16.0
2 HW 1 SOL HW1 1 0.417 1.008
3 HW 1 SOL HW2 1 0.417 1.008
[ bonds ]
1 2 1
1 3 1
[ system ]
case
[ molecules ]
PEP 1
SOL 2
"""
    )

    system, mols = Xponge.load_molitp(str(top), water_replace=False)
    current_summary = {
        "system_name": system.name,
        "atom_count": system.atom_count,
        "residue_names": [res.name for res in system.residues],
        "residue_counts": system.residue_counts(),
        "mols": {name: [res.name for res in mol.residues] for name, mol in mols.items()},
    }

    script = textwrap.dedent(
        f"""
        import json
        import sys
        sys.path.insert(0, {str(XPONGE_REPO)!r})
        import Xponge
        import Xponge.forcefield.amber.ff14sb

        system, mols = Xponge.load_molitp({str(top)!r}, water_replace=False)
        payload = {{
            "system_name": system.name,
            "atom_count": len(system.atoms),
            "residue_names": [res.name for res in system.residues],
            "residue_counts": {{}},
            "mols": {{}},
        }}
        for residue in system.residues:
            payload["residue_counts"][residue.name] = payload["residue_counts"].get(residue.name, 0) + 1
        for name, mol in mols.items():
            if hasattr(mol, "residues"):
                payload["mols"][name] = [res.name for res in mol.residues]
            else:
                payload["mols"][name] = [mol.name]
        print(json.dumps(payload, sort_keys=True))
        """
    )
    reference = json.loads(_run_original_xponge(script))

    assert current_summary == reference


def test_load_molitp_water_replace_maps_sol_to_wat_when_tip3p_is_loaded(tmp_path):
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    top = tmp_path / "system.top"
    top.write_text(
        """[ moleculetype ]
SOL 3
[ atoms ]
1 OW 1 SOL OW 1 -0.834 16.0
2 HW 1 SOL HW1 1 0.417 1.008
3 HW 1 SOL HW2 1 0.417 1.008
[ bonds ]
1 2 1
1 3 1
[ system ]
water
[ molecules ]
SOL 2
"""
    )

    system, mols = Xponge.load_molitp(str(top), water_replace=True)

    assert sorted(mols) == ["SOL"]
    assert [res.name for res in mols["SOL"].residues] == ["WAT"]
    assert system.residue_counts() == {"WAT": 2}


def test_load_molitp_invokes_registered_bonded_type_parsers_and_copies_special_forces(tmp_path):
    top = tmp_path / "system.top"
    top.write_text(
        """[ moleculetype ]
MOL 3
[ atoms ]
1 Q 1 LIG A1 1 0.0 12.0
2 Q 1 LIG A2 1 0.0 12.0
3 Q 1 LIG A3 1 0.0 12.0
[ bonds ]
1 2 1
2 3 1
[ pairs ]
1 3 1
[ system ]
case
[ molecules ]
MOL 2
"""
    )

    def pair_parser(words, mol, _stat):
        mol.add_nb14_extra(int(words[0]) - 1, int(words[1]) - 1, 1.25, 2.5, 0.75)

    Xponge.register_amber_lj_parameter("Q", "Q", 0.2, 1.0)
    Xponge.GlobalSetting.Set_GMX_Bonded_Type_Parser("pair", 1, pair_parser)
    try:
        system, mols = Xponge.load_molitp(str(top), water_replace=False)
    finally:
        Xponge.GlobalSetting._gmx_bonded_type_parsers.pop(("pair", 1), None)

    assert sorted(mols) == ["MOL"]
    out = Xponge.Save_SPONGE_Input(system, prefix="molitp", dirname=str(tmp_path))
    assert "nb14_extra" in out
    lines = (tmp_path / "molitp_nb14_extra.txt").read_text().splitlines()
    assert lines[0] == "2"
    records = [[float(value) for value in line.split()] for line in lines[1:]]
    assert records == [
        [0.0, 2.0, 1.25, 2.5, 0.75],
        [3.0, 5.0, 1.25, 2.5, 0.75],
    ]


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


def test_opls_itp_parser_supports_larger_loader_driven_exports_with_pairs_and_rb(tmp_path):
    itp = tmp_path / "opls_big.itp"
    itp.write_text(
        """[ atomtypes ]
opls_1 C 6 12.011 0.0 A 0.35 0.276144
opls_2 H 1 1.008 0.0 A 0.25 0.125520
[ moleculetype ]
MOL 3
[ atoms ]
1 opls_1 1 MOL A1 1 0.0 12.011
2 opls_1 1 MOL A2 1 0.0 12.011
3 opls_1 1 MOL A3 1 0.0 12.011
4 opls_1 1 MOL A4 1 0.0 12.011
5 opls_2 1 MOL H5 1 0.0 1.008
[ bonds ]
1 2 1
2 3 1
3 4 1
4 5 1
[ pairs ]
1 4 1 1.25 2.5 0.75
2 5 1 1.15 2.2 0.55
[ dihedrals ]
1 2 3 4 3 0.1 0.2 0.3 0.4 0.5 0.6
2 3 4 5 3 0.6 0.5 0.4 0.3 0.2 0.1
"""
    )

    mol = Xponge.load_opls_itp_file(str(itp))
    out = Xponge.Save_SPONGE_Input(mol, prefix="opls_big", dirname=str(tmp_path))

    assert mol.atom_count == 5
    assert mol.residue_count == 1
    assert mol.validate()
    assert {"Ryckaert_Bellemans", "nb14_extra", "listed_forces"}.issubset(out)
    assert (tmp_path / "opls_big_bond.txt").read_text().splitlines()[0] == "4"
    assert (tmp_path / "opls_big_angle.txt").read_text().splitlines()[0] == "3"
    assert (tmp_path / "opls_big_nb14.txt").read_text().splitlines()[0] == "2"
    assert (tmp_path / "opls_big_nb14_extra.txt").read_text().splitlines() == [
        "2",
        "0 3 1.250000e+00 2.500000e+00 7.500000e-01",
        "1 4 1.150000e+00 2.200000e+00 5.500000e-01",
    ]
    assert (tmp_path / "opls_big_Ryckaert_Bellemans.txt").read_text().splitlines() == [
        "2",
        "0 1 2 3 0.100000 0.200000 0.300000 0.400000 0.500000 0.600000",
        "1 2 3 4 0.600000 0.500000 0.400000 0.300000 0.200000 0.100000",
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


def test_charmm_prm_and_rtf_parser_supports_larger_loader_driven_export(tmp_path):
    prm = tmp_path / "larger.prm"
    prm.write_text(
        """MASS -1 XCT 12.011
MASS -1 XHT 1.008
NONBONDED
XCT 0.0 -0.100 2.000
XHT 0.0 -0.050 1.000
ANGLE
XHT XCT XCT 45.0 112.0 5.0 2.2
XCT XCT XHT 45.0 112.0 5.0 2.2
XHT XCT XHT 50.0 109.5 10.0 2.0
NBFIX
XCT XHT 1.25 2.5 0.75
"""
    )
    rtf = tmp_path / "larger.rtf"
    rtf.write_text(
        """RESI BUT 0.0
ATOM H1 XHT 0.1
ATOM C1 XCT -0.2
ATOM C2 XCT -0.2
ATOM H2 XHT 0.1
BOND H1 C1 C1 C2 C2 H2
END
"""
    )

    Xponge.load_charmm_parameter_file(str(prm))
    mol = Xponge.load_charmm_topology_file(str(rtf))
    out = Xponge.Save_SPONGE_Input(mol, prefix="charmm_big", dirname=str(tmp_path))

    assert mol.residue_count == 1
    assert mol.atom_count == 4
    assert [atom.name for atom in mol.residues[0].atoms] == ["H1", "C1", "C2", "H2"]
    assert mol.validate()
    assert {"urey_bradley", "nb14_extra"}.issubset(out)
    assert (tmp_path / "charmm_big_resname.txt").read_text().splitlines() == ["1", "BUT"]
    assert (tmp_path / "charmm_big_bond.txt").read_text().splitlines()[0] == "3"
    assert (tmp_path / "charmm_big_angle.txt").read_text().splitlines()[0] == "2"
    assert (tmp_path / "charmm_big_nb14.txt").read_text().splitlines()[0] == "1"
    assert (tmp_path / "charmm_big_urey_bradley.txt").read_text().splitlines() == [
        "2",
        "0 1 2 45.000000 112.000000 5.000000 2.200000",
        "1 2 3 45.000000 112.000000 5.000000 2.200000",
    ]
    assert (tmp_path / "charmm_big_nb14_extra.txt").read_text().splitlines() == [
        "4",
        "0 1 1.250000e+00 2.500000e+00 0.000000e+00",
        "0 2 1.250000e+00 2.500000e+00 0.000000e+00",
        "1 3 1.250000e+00 2.500000e+00 0.000000e+00",
        "2 3 1.250000e+00 2.500000e+00 0.000000e+00",
    ]


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


def test_martini300_module_supports_loader_driven_small_molecule_export_workflow(tmp_path):
    import XpongeCPP.forcefield.martini.martini300  # noqa: F401

    top = tmp_path / "chex.top"
    top.write_text(
        """#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_small_molecules_v1.itp"
[ system ]
martini_chex
[ molecules ]
CHEX 1
"""
    )

    system, mols = Xponge.load_molitp(str(top), water_replace=False)
    out = Xponge.Save_SPONGE_Input(system, prefix="martini_chex", dirname=str(tmp_path))

    assert system.name == "martini_chex"
    assert [res.name for res in system.residues] == ["CHEX"]
    assert sorted(mols) == ["CHEX"]
    assert system.atom_count == 2
    assert system.validate()
    assert (tmp_path / "martini_chex_resname.txt").read_text().splitlines() == ["1", "CHEX"]
    assert (tmp_path / "martini_chex_bond.txt").read_text().splitlines()[0] == "1"
    assert "LJ" in out


def test_martini300_constraints_topology_reports_current_connectivity_limitation_explicitly(tmp_path):
    import XpongeCPP.forcefield.martini.martini300  # noqa: F401

    top = tmp_path / "4mimi.top"
    top.write_text(
        """#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_small_molecules_v1.itp"
[ system ]
martini_4mimi
[ molecules ]
4MIMI 1
"""
    )

    system, mols = Xponge.load_molitp(str(top), water_replace=False)

    assert system.atom_count == 3
    assert [res.name for res in system.residues] == ["4MIMI"]
    assert sorted(mols) == ["4MIMI"]
    with pytest.raises(RuntimeError, match="missing residue template/connectivity for residue: 4MIMI"):
        Xponge.Save_SPONGE_Input(system, prefix="martini_4mimi", dirname=str(tmp_path))


def test_non_amber_forcefield_import_modules_are_available():
    import XpongeCPP.forcefield.charmm.charmm36  # noqa: F401
    import XpongeCPP.forcefield.martini.martini300  # noqa: F401
    import XpongeCPP.forcefield.opls.oplsaam  # noqa: F401
    import XpongeCPP.forcefield.sw.mw as mw
    import XpongeCPP.forcefield.edip.si as si

    assert callable(mw.load_parameters)
    assert callable(si.load_parameters)


def test_charmm36_module_registers_protein_dna_and_rna_templates_with_reference_terminal_mappings():
    import XpongeCPP.forcefield.charmm.charmm36  # noqa: F401

    protein_pdb = "".join([
        _pdb_atom_line(1, "N", "HIS", "A", 1, 0.0, 0.0, 0.0, "N"),
        _pdb_atom_line(2, "CA", "HIS", "A", 1, 1.4, 0.0, 0.0, "C"),
        _pdb_atom_line(3, "C", "HIS", "A", 1, 2.8, 0.0, 0.0, "C"),
        _pdb_atom_line(4, "O", "HIS", "A", 1, 3.9, 0.0, 0.0, "O"),
        _pdb_atom_line(5, "HD1", "HIS", "A", 1, 1.2, 1.1, 0.0, "H"),
        _pdb_atom_line(6, "N", "CYS", "A", 2, 4.2, 0.0, 0.0, "N"),
        _pdb_atom_line(7, "CA", "CYS", "A", 2, 5.6, 0.0, 0.0, "C"),
        _pdb_atom_line(8, "C", "CYS", "A", 2, 7.0, 0.0, 0.0, "C"),
        _pdb_atom_line(9, "O", "CYS", "A", 2, 8.1, 0.0, 0.0, "O"),
        _pdb_atom_line(10, "OXT", "CYS", "A", 2, 7.0, 1.2, 0.0, "O"),
        "TER\n",
        "END\n",
    ])
    dna_pdb = "".join([
        _pdb_atom_line(1, "O5'", "DA", "A", 1, -0.6, -1.0, 0.0, "O"),
        _pdb_atom_line(2, "C5'", "DA", "A", 1, -0.1, -2.2, 0.0, "C"),
        _pdb_atom_line(3, "C3'", "DA", "A", 1, -1.8, -0.3, 0.0, "C"),
        _pdb_atom_line(4, "O3'", "DA", "A", 1, -2.9, -0.9, 0.0, "O"),
        _pdb_atom_line(5, "P", "DT", "A", 2, -3.9, -0.1, 0.0, "P"),
        _pdb_atom_line(6, "OP2", "DT", "A", 2, -5.3, -0.1, 0.0, "O"),
        _pdb_atom_line(7, "O5'", "DT", "A", 2, -3.4, 1.2, 0.0, "O"),
        _pdb_atom_line(8, "C5'", "DT", "A", 2, -3.9, 2.4, 0.0, "C"),
        _pdb_atom_line(9, "C3'", "DT", "A", 2, -2.2, 1.8, 0.0, "C"),
        _pdb_atom_line(10, "O3'", "DT", "A", 2, -1.2, 1.1, 0.0, "O"),
        "TER\n",
        "END\n",
    ])
    rna_pdb = "".join([
        _pdb_atom_line(1, "O5'", "A", "A", 1, -0.6, 4.0, 0.0, "O"),
        _pdb_atom_line(2, "C5'", "A", "A", 1, -0.1, 2.8, 0.0, "C"),
        _pdb_atom_line(3, "C3'", "A", "A", 1, -1.8, 4.7, 0.0, "C"),
        _pdb_atom_line(4, "O3'", "A", "A", 1, -2.9, 4.1, 0.0, "O"),
        _pdb_atom_line(5, "P", "U", "A", 2, -3.9, 4.9, 0.0, "P"),
        _pdb_atom_line(6, "OP2", "U", "A", 2, -5.3, 4.9, 0.0, "O"),
        _pdb_atom_line(7, "O5'", "U", "A", 2, -3.4, 6.2, 0.0, "O"),
        _pdb_atom_line(8, "C5'", "U", "A", 2, -3.9, 7.4, 0.0, "C"),
        _pdb_atom_line(9, "C3'", "U", "A", 2, -2.2, 6.8, 0.0, "C"),
        _pdb_atom_line(10, "O3'", "U", "A", 2, -1.2, 6.1, 0.0, "O"),
        "TER\n",
        "END\n",
    ])

    current = {
        "templates": {
            "NALA": Xponge.template_atom_count("NALA"),
            "CALA": Xponge.template_atom_count("CALA"),
            "HIS": Xponge.template_atom_count("HIS"),
            "DA5": Xponge.template_atom_count("DA5"),
            "DT3": Xponge.template_atom_count("DT3"),
            "A5": Xponge.template_atom_count("A5"),
            "U3": Xponge.template_atom_count("U3"),
        },
        "protein_residues": [res.name for res in Xponge.load_pdb(StringIO(protein_pdb)).residues],
        "dna_residues": [res.name for res in Xponge.load_pdb(StringIO(dna_pdb)).residues],
        "rna_residues": [res.name for res in Xponge.load_pdb(StringIO(rna_pdb)).residues],
    }

    script = textwrap.dedent(
        f"""
        import json
        import sys
        from io import StringIO
        sys.path.insert(0, {str(XPONGE_REPO)!r})
        import Xponge
        import Xponge.forcefield.charmm.charmm36

        protein_pdb = {protein_pdb!r}
        dna_pdb = {dna_pdb!r}
        rna_pdb = {rna_pdb!r}
        payload = {{
            "templates": {{}},
            "protein_residues": [res.name for res in Xponge.load_pdb(StringIO(protein_pdb)).residues],
            "dna_residues": [res.name for res in Xponge.load_pdb(StringIO(dna_pdb)).residues],
            "rna_residues": [res.name for res in Xponge.load_pdb(StringIO(rna_pdb)).residues],
        }}
        for name in ("NALA", "CALA", "HIS", "DA5", "DT3", "A5", "U3"):
            payload["templates"][name] = len(Xponge.ResidueType.get_type(name).atoms)
        print(json.dumps(payload, sort_keys=True))
        """
    )
    reference = json.loads(_run_original_xponge(script))

    assert current == reference


def test_oplsaam_module_registers_protein_templates_with_reference_terminal_mappings():
    import XpongeCPP.forcefield.opls.oplsaam  # noqa: F401

    protein_pdb = "".join([
        _pdb_atom_line(1, "N", "HIS", "A", 1, 0.0, 0.0, 0.0, "N"),
        _pdb_atom_line(2, "CA", "HIS", "A", 1, 1.4, 0.0, 0.0, "C"),
        _pdb_atom_line(3, "C", "HIS", "A", 1, 2.8, 0.0, 0.0, "C"),
        _pdb_atom_line(4, "O", "HIS", "A", 1, 3.9, 0.0, 0.0, "O"),
        _pdb_atom_line(5, "HD1", "HIS", "A", 1, 1.2, 1.1, 0.0, "H"),
        _pdb_atom_line(6, "N", "CYS", "A", 2, 4.2, 0.0, 0.0, "N"),
        _pdb_atom_line(7, "CA", "CYS", "A", 2, 5.6, 0.0, 0.0, "C"),
        _pdb_atom_line(8, "C", "CYS", "A", 2, 7.0, 0.0, 0.0, "C"),
        _pdb_atom_line(9, "O", "CYS", "A", 2, 8.1, 0.0, 0.0, "O"),
        _pdb_atom_line(10, "OXT", "CYS", "A", 2, 7.0, 1.2, 0.0, "O"),
        "TER\n",
        "END\n",
    ])

    current = {
        "templates": {
            "NALA": Xponge.template_atom_count("NALA"),
            "CALA": Xponge.template_atom_count("CALA"),
            "HIS": Xponge.template_atom_count("HIS"),
        },
        "protein_residues": [res.name for res in Xponge.load_pdb(StringIO(protein_pdb)).residues],
    }

    script = textwrap.dedent(
        f"""
        import json
        import sys
        from io import StringIO
        sys.path.insert(0, {str(XPONGE_REPO)!r})
        import Xponge
        import Xponge.forcefield.opls.oplsaam

        protein_pdb = {protein_pdb!r}
        payload = {{
            "templates": {{}},
            "protein_residues": [res.name for res in Xponge.load_pdb(StringIO(protein_pdb)).residues],
        }}
        for name in ("NALA", "CALA", "HIS"):
            payload["templates"][name] = len(Xponge.ResidueType.get_type(name).atoms)
        print(json.dumps(payload, sort_keys=True))
        """
    )
    reference = json.loads(_run_original_xponge(script))

    assert current == reference


def test_charmm36_template_peptide_exports_representative_common_files(tmp_path):
    import XpongeCPP.forcefield.charmm.charmm36  # noqa: F401

    mol = Xponge.get_template_molecule("NALA") + Xponge.get_template_molecule("CGLY")

    assert [res.name for res in mol.residues] == ["NALA", "CGLY"]
    assert mol.validate()

    out = Xponge.Save_SPONGE_Input(mol, prefix="charmm36_peptide", dirname=str(tmp_path))

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
    assert (tmp_path / "charmm36_peptide_bond.txt").read_text().splitlines()[0] == "19"
    assert (tmp_path / "charmm36_peptide_resname.txt").read_text().splitlines() == ["2", "NALA", "CGLY"]


def test_oplsaam_template_peptide_exports_representative_common_files(tmp_path):
    import XpongeCPP.forcefield.opls.oplsaam  # noqa: F401

    mol = Xponge.get_template_molecule("NALA") + Xponge.get_template_molecule("CGLY")

    assert [res.name for res in mol.residues] == ["NALA", "CGLY"]
    assert mol.validate()

    out = Xponge.Save_SPONGE_Input(mol, prefix="oplsaam_peptide", dirname=str(tmp_path))

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
    assert (tmp_path / "oplsaam_peptide_bond.txt").read_text().splitlines()[0] == "19"
    assert (tmp_path / "oplsaam_peptide_resname.txt").read_text().splitlines() == ["2", "NALA", "CGLY"]


def test_martini300_module_supports_loader_driven_solvent_export_workflow(tmp_path):
    import XpongeCPP.forcefield.martini.martini300  # noqa: F401

    top = tmp_path / "martini.top"
    top.write_text(
        """#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
[ system ]
martini_dmso
[ molecules ]
DMSO 1
"""
    )

    system, mols = Xponge.load_molitp(str(top), water_replace=False)

    assert system.name == "martini_dmso"
    assert sorted(mols) == ["DMSO"]
    assert [res.name for res in system.residues] == ["DMS"]
    assert system.validate()

    out = Xponge.Save_SPONGE_Input(system, prefix="martini_dmso", dirname=str(tmp_path))

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
    assert (tmp_path / "martini_dmso_bond.txt").read_text().splitlines()[0] == "1"
    assert (tmp_path / "martini_dmso_resname.txt").read_text().splitlines() == ["1", "DMS"]
