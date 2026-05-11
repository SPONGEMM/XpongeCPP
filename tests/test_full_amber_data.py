from importlib import resources
from io import StringIO

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


def test_reference_xponge_forcefield_tree_is_packaged():
    reference_root = resources.files("XpongeCPP").joinpath("data", "reference_forcefield")

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
