from importlib import resources

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
