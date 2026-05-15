import XpongeCPP as Xponge
import pytest


def test_compat_module_adds_instance_style_molecule_save_methods(tmp_path):
    import XpongeCPP.compat  # noqa: F401
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.get_template_molecule("ALA")

    pdb_path = tmp_path / "lig.pdb"
    mol2_path = tmp_path / "lig.mol2"
    sponge_dir = tmp_path / "spg"
    sponge_dir.mkdir()

    mol.save_pdb(str(pdb_path))
    mol.save_mol2(str(mol2_path))
    outputs = mol.save_sponge_input(prefix="lig", dirname=str(sponge_dir))

    assert pdb_path.is_file()
    assert mol2_path.is_file()
    assert outputs is mol
    assert (sponge_dir / "lig_bond.txt").is_file()
    assert (sponge_dir / "lig_bond.txt").is_file()


def test_compat_layer_can_inject_legacy_template_globals():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    from XpongeCPP.compat import enable_legacy_namespace

    namespace = {}
    enable_legacy_namespace(namespace, template_names=["NALA", "ALA", "CALA"])

    mol = namespace["NALA"] + namespace["ALA"] * 12 + namespace["CALA"]

    assert mol.residue_count == 14
    assert [res.name for res in mol.residues[:2]] == ["NALA", "ALA"]
    assert mol.residues[-1].name == "CALA"
    assert mol.validate()


def test_legacy_residuetype_handles_support_add_and_mul():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = (
        Xponge.ResidueType.get_type("NALA")
        + Xponge.ResidueType.get_type("ALA") * 3
        + Xponge.ResidueType.get_type("CALA")
    )

    assert mol.residue_count == 5
    assert [res.name for res in mol.residues] == ["NALA", "ALA", "ALA", "ALA", "CALA"]
    assert mol.validate()


def test_legacy_ff19sb_gb_workflow_shim_runs(tmp_path):
    import XpongeCPP.forcefield.amber.ff19sb  # noqa: F401
    from XpongeCPP.forcefield.special import gb

    mol = (
        Xponge.ResidueType.get_type("NALA")
        + Xponge.ResidueType.get_type("ALA") * 10
        + Xponge.ResidueType.get_type("CALA")
    )
    gb.set_gb_radius(mol)

    Xponge.Save_PDB(mol, tmp_path / "ALA.pdb")
    outputs = Xponge.Save_SPONGE_Input(mol, prefix="ALA", dirname=str(tmp_path))

    assert mol.residue_count == 12
    assert mol.atom_count > 0
    assert (tmp_path / "ALA.pdb").is_file()
    assert (tmp_path / "ALA_gb.txt").is_file()
    assert outputs is mol


def test_xponge_package_alias_supports_ff19sb_gb_workflow(tmp_path, monkeypatch):
    import sys
    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401
    from Xponge.forcefield.special import gb

    monkeypatch.chdir(tmp_path)
    sys.modules.pop("Xponge", None)

    mol = (
        Xponge.ResidueType.get_type("NALA")
        + Xponge.ResidueType.get_type("ALA") * 2
        + Xponge.ResidueType.get_type("CALA")
    )
    gb.set_gb_radius(mol)

    Xponge.Save_PDB(mol, "alias.pdb")
    outputs = Xponge.Save_SPONGE_Input(mol, "alias")

    assert mol.residue_count == 4
    assert (tmp_path / "alias.pdb").is_file()
    assert (tmp_path / "alias_gb.txt").is_file()
    assert outputs is mol


def test_xpongecpp_common_legacy_import_paths_resolve():
    import XpongeCPP.assign as assign
    import XpongeCPP.build as build
    import XpongeCPP.forcefield.base.angle_base as angle_base
    import XpongeCPP.forcefield.base.bond_base as bond_base
    import XpongeCPP.forcefield.base.charge_base as charge_base
    import XpongeCPP.forcefield.base.exclude_base as exclude_base
    import XpongeCPP.forcefield.base.lj_base as lj_base
    import XpongeCPP.forcefield.base.mass_base as mass_base
    import XpongeCPP.load as load
    from XpongeCPP.helper import GlobalSetting, Xdict, Xprint

    assert load.GromacsTopologyIterator is Xponge.GromacsTopologyIterator
    assert callable(build.save_pdb)
    assert assign.AssignRule is Xponge.AssignRule
    assert issubclass(Xdict, dict)
    assert callable(Xprint)
    assert GlobalSetting is Xponge.GlobalSetting
    assert angle_base.AngleType is not None
    assert bond_base.BondType is not None
    assert lj_base.LJType is not None
    assert exclude_base.Exclude(4).n == 4
    assert charge_base is not None
    assert mass_base is not None


def test_xponge_package_alias_supports_common_legacy_import_paths():
    import Xponge.assign as assign
    import Xponge.build as build
    import Xponge.forcefield.base.angle_base as angle_base
    import Xponge.forcefield.base.bond_base as bond_base
    import Xponge.forcefield.base.charge_base as charge_base
    import Xponge.forcefield.base.exclude_base as exclude_base
    import Xponge.forcefield.base.lj_base as lj_base
    import Xponge.forcefield.base.mass_base as mass_base
    import Xponge.load as load
    from Xponge.helper import GlobalSetting, Xdict, Xpri

    assert load.GromacsTopologyIterator is Xponge.GromacsTopologyIterator
    assert callable(build.save_pdb)
    assert assign.AssignRule is Xponge.AssignRule
    assert issubclass(Xdict, dict)
    assert callable(Xpri)
    assert GlobalSetting is Xponge.GlobalSetting
    assert angle_base.AngleType is not None
    assert bond_base.BondType is not None
    assert lj_base.LJType is not None
    assert exclude_base.Exclude(4).n == 4
    assert charge_base is not None
    assert mass_base is not None


def test_xponge_package_alias_supports_high_frequency_forcefield_and_helper_modules():
    import Xponge.forcefield.amber.bsc1 as bsc1
    import Xponge.forcefield.amber.ol3 as ol3
    import Xponge.forcefield.amber.tip4pew as tip4pew
    import Xponge.forcefield.charmm.charmm27 as charmm27
    import Xponge.forcefield.charmm.charmm36 as charmm36
    import Xponge.forcefield.charmm.tip3p as charmm_tip3p
    import Xponge.forcefield.charmm.tip3p_charmm as tip3p_charmm
    import Xponge.forcefield.opls.oplsaam as oplsaam
    import Xponge.forcefield.special.min as min_helper
    import Xponge.forcefield.sw.mw as mw
    from Xponge.helper.gromacs import Sort_Atoms_By_Gro

    assert bsc1 is not None
    assert ol3 is not None
    assert tip4pew is not None
    assert charmm27 is not None
    assert charmm36 is not None
    assert charmm_tip3p is not None
    assert tip3p_charmm is not None
    assert oplsaam is not None
    assert min_helper is not None
    assert mw is not None
    assert callable(Sort_Atoms_By_Gro)


def test_xponge_package_alias_special_min_helpers_toggle_fake_parameter_exports(tmp_path, monkeypatch):
    from Xponge.forcefield.special import min as min_helper
    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401

    assert hasattr(min_helper, "save_min_bonded_parameters")
    assert hasattr(min_helper, "do_not_save_min_bonded_parameters")
    assert hasattr(min_helper, "Save_Min_Bonded_Parameters")
    assert hasattr(min_helper, "Do_Not_Save_Min_Bonded_Parameters")

    monkeypatch.chdir(tmp_path)
    mol = Xponge.NALA + Xponge.ALA + Xponge.CALA

    min_helper.Save_Min_Bonded_Parameters()
    outputs = Xponge.Save_SPONGE_Input(mol, "mincase")
    assert outputs is mol
    assert (tmp_path / "mincase_fake_mass.txt").is_file()
    assert (tmp_path / "mincase_fake_LJ.txt").is_file()
    assert (tmp_path / "mincase_fake_charge.txt").is_file()

    min_helper.Do_Not_Save_Min_Bonded_Parameters()
    Xponge.Save_SPONGE_Input(mol, "nomincase")
    assert not (tmp_path / "nomincase_fake_mass.txt").exists()
    assert not (tmp_path / "nomincase_fake_LJ.txt").exists()
    assert not (tmp_path / "nomincase_fake_charge.txt").exists()


def test_xponge_package_alias_special_fep_module_exports_legacy_surface():
    from Xponge.forcefield.special import fep

    names = [
        "prepare_lj_soft_core",
        "Prepare_LJ_Soft_Core",
        "set_lj_type_b",
        "Set_LJ_Type_B",
        "set_subsys",
        "Set_Subsys",
        "enable_lj_soft_core",
        "Enable_LJ_Soft_Core",
        "merge_dual_topology",
        "Merge_Dual_Topology",
        "merge_force_field",
        "Merge_Force_Field",
        "get_free_molecule",
        "Get_Free_Molecule",
        "intramolecule_nb_to_nb14",
        "Intramolecule_NB_To_NB14",
        "save_soft_core_lj",
        "Save_Soft_Core_LJ",
        "save_hard_core_lj",
        "Save_Hard_Core_LJ",
        "add_soft_bond_from_a",
        "Add_Soft_Bond_From_A",
        "add_soft_bond_from_b",
        "Add_Soft_Bond_From_B",
    ]
    for name in names:
        assert hasattr(fep, name), name


def test_xponge_package_alias_supports_minimal_cvsystem_workflow(tmp_path):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    from Xponge.helper.cv import CVSystem

    mol = (
        Xponge.ResidueType.get_type("ACE")
        + Xponge.ResidueType.get_type("ALA") * 2
        + Xponge.ResidueType.get_type("NME")
    )
    cv = CVSystem(mol)
    cv.add_center("c", "protein")
    cv.add_cv_position("x", "c", "x", scaled=False)
    cv.print("x")
    cv.add_cv_rmsd("r", "protein and backbone")
    cv.print("r")
    cv.steer("x", 2)
    cv.output(tmp_path / "cv.txt", folder=tmp_path)

    assert (tmp_path / "cv.txt").is_file()
    assert (tmp_path / "r_atom.txt").is_file()
    assert (tmp_path / "r_coordinate.txt").is_file()
    text = (tmp_path / "cv.txt").read_text()
    assert "definition of virtual atoms" in text
    assert "CV_type = position_x" in text
    assert "CV_type = rmsd" in text
    assert "steer" in text


def test_xponge_package_alias_supports_cvsystem_meta1d_output(tmp_path):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    from Xponge.helper.cv import CVSystem

    mol = Xponge.ResidueType.get_type("ACE") + Xponge.ResidueType.get_type("ALA") + Xponge.ResidueType.get_type("NME")
    cv = CVSystem(mol)
    cv.add_cv_dihedral("torsion", mol.atoms[2], mol.atoms[0], mol.atoms[1], mol.atoms[3])
    cv.meta1d("torsion", CV_grid=100, CV_minimal=-3.14, CV_maximum=3.14, CV_sigma=0.5, welltemp_factor=50, height=1)
    cv.print("torsion")
    cv.output(tmp_path / "cv_meta.txt", folder=tmp_path)

    text = (tmp_path / "cv_meta.txt").read_text()
    assert "metad" in text
    assert "CV = torsion" in text
    assert "CV_grid = 100" in text
    assert "CV_period" in text


def test_xponge_top_level_and_package_style_imports_match_legacy_shape():
    from Xponge import (
        AbstractMolecule,
        Assign,
        Atom,
        AtomType,
        Entity,
        GlobalSetting,
        ResidueLink,
        Type,
        Xpri,
        Xopen,
        bar,
        kb,
        pdb_filter,
        pi,
    )
    from Xponge.analysis import MdoutReader, wham
    from Xponge.analysis.md_analysis import XpongeMoleculeReader, mda
    from Xponge.analysis.sasa import SASA
    from Xponge.forcefield.amber import gaff
    from Xponge.helper.file import file_filter
    from Xponge.helper.math import get_fibonacci_grid
    from Xponge.helper.namespace import source as namespace_source
    from Xponge.load import GromacsTopologyIterator
    from Xponge.mdrun import run
    from Xponge.tools.unittests import CATEGORY, XpongeTestRunner, mytest

    assert GlobalSetting is not None
    assert callable(Xpri)
    assert callable(Xopen)
    assert callable(pdb_filter)
    assert callable(file_filter)
    assert callable(get_fibonacci_grid)
    assert callable(namespace_source)
    assert pi is not None
    assert kb is not None
    assert bar is not None
    assert AbstractMolecule is not None
    assert Atom is not None
    assert Assign is not None
    assert AtomType is not None
    assert Entity is not None
    assert ResidueLink is not None
    assert Type is not None
    assert gaff is not None
    assert GromacsTopologyIterator is not None
    assert MdoutReader is not None
    assert hasattr(wham, "WHAM")
    assert XpongeMoleculeReader is not None
    assert mda is None
    assert SASA is not None
    assert callable(run)
    assert CATEGORY["0"] == "base"
    assert XpongeTestRunner is not None
    assert callable(mytest)


def test_legacy_analysis_and_mdrun_surfaces_have_first_wave_real_and_placeholder_behavior(tmp_path):
    import numpy as np
    from Xponge.analysis import MdoutReader, wham
    from Xponge.analysis.md_analysis import XpongeMoleculeReader, mda
    from Xponge.analysis.sasa import SASA
    from Xponge.mdrun import run

    mdout_path = tmp_path / "demo.mdout"
    mdout_path.write_text(
        "step Time Temperature Potential LJ\n"
        "0 0 300.0 100.0 10.0\n"
        "1 0.1 299.8 101.1 11.1\n",
        encoding="utf-8",
    )
    mdout = MdoutReader(mdout_path)
    assert list(mdout.step) == [0.0, 1.0]
    assert list(mdout.Potential) == [100.0, 101.1]

    w = wham.WHAM(np.linspace(-1.0, 1.0, 5), 300.0, 10.0, np.linspace(-0.5, 0.5, 4))
    assert w.window_edges.shape == (5,)

    assert mda is None
    for callable_obj in (XpongeMoleculeReader, SASA):
        try:
            callable_obj()
        except (NotImplementedError, ModuleNotFoundError) as exc:
            message = str(exc).lower()
            assert "compatibility" in message or "legacy" in message or "analysis" in message or "mdanalysis" in message
        else:  # pragma: no cover
            raise AssertionError(f"{callable_obj} should fail with a clear NotImplementedError")

    try:
        run(["mdrun", "-reset"])
    except SystemExit as exc:
        assert exc.code in (None, 0)
    else:  # pragma: no cover
        raise AssertionError("mdrun -reset should exit after updating BIN_PATH.dat")


def test_old_xponge_init_export_surface_has_first_wave_compatibility_shape():
    import Xponge
    from Xponge import build_bonded_force, get_mindsponge_system_energy, set_global_alternative_names

    assert callable(set_global_alternative_names)
    assert callable(build_bonded_force)
    assert callable(get_mindsponge_system_energy)
    assert hasattr(Xponge, "Atom")
    assert hasattr(Xponge, "Residue")


def test_mindsponge_registry_surface_is_dependency_aware_and_no_longer_placeholder():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    assert hasattr(Xponge.Molecule, "Set_MindSponge_Todo")
    assert hasattr(Xponge.Molecule, "Del_MindSponge_Todo")
    assert hasattr(Xponge.Molecule, "_mindsponge_todo")
    assert hasattr(Xponge.Molecule, "Set_Save_SPONGE_Input")
    assert hasattr(Xponge.Molecule, "_save_functions")
    assert "coordinate" in Xponge.Molecule._mindsponge_todo
    assert "mass" in Xponge.Molecule._mindsponge_todo
    assert "charge" in Xponge.Molecule._mindsponge_todo

    mol = Xponge.ACE + Xponge.ALA + Xponge.NME
    try:
        result = Xponge.get_mindsponge_system_energy(mol)
    except ModuleNotFoundError as exc:
        assert "mindsponge" in str(exc).lower() or "mindspore" in str(exc).lower()
    else:  # pragma: no cover - optional dependency path
        system, energy = result
        assert system is not None
        assert energy is not None


def test_build_bonded_force_is_first_wave_runnable_compatibility_for_core_objects(tmp_path):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    import Xponge.build as build_mod

    mol = Xponge.ACE + Xponge.ALA + Xponge.NME
    residue = mol.residues[1]
    residue_type = Xponge.ResidueType.get_type("ALA")

    assert Xponge.build_bonded_force(mol) is mol
    assert Xponge.build_bonded_force(residue) is residue
    assert Xponge.build_bonded_force(residue_type) is residue_type
    assert build_mod.build_bonded_force(mol) is mol

    outputs = Xponge.Save_SPONGE_Input(Xponge.build_bonded_force(mol), prefix="built", dirname=str(tmp_path))
    assert outputs is mol
    assert (tmp_path / "built_bond.txt").is_file()


def test_legacy_atomtype_registry_supports_old_style_new_from_string_and_get_type():
    import Xponge

    Xponge.AtomType.New_From_String(
        """
name  mass   charge[e]
H     1.008  +0.25
C     12.00  -1.00
"""
    )

    h = Xponge.AtomType.get_type("H")
    c = Xponge.AtomType.Get_Type("C")
    all_types = Xponge.AtomType.get_all_types()

    assert h.name == "H"
    assert h.mass == 1.008
    assert h.__dict__["charge[e]"] == 0.25
    assert c.name == "C"
    assert c.mass == 12.0
    assert c.__dict__["charge[e]"] == -1.0
    assert all_types["H"] is h
    assert all_types["C"] is c


def test_legacy_basic_manual_building_shape_matches_old_test0_base_usage():
    import Xponge

    assign = Xponge.Assign()
    assign.add_atom("O", 0, 0, 0)
    assign.addAtom("H1", 0, 1, 0)
    assign.Add_Atom("H2", 1, 0, 0)
    assign.add_bond(0, 1, 1)
    assign.addBond(0, 2, 1)
    assign.Add_Bond(1, 2, 1)

    Xponge.AtomType.New_From_String(
        """
name  mass   charge[e]
H     1.008  +0.25
C     12.00  -1.00
"""
    )
    h = Xponge.AtomType.get_type("H")
    c = Xponge.AtomType.get_type("C")
    xyj = Xponge.ResidueType(name="XYJ")
    xyj.add_atom("H1", h, 1, 0, 0)
    xyj.add_atom("H2", "H", 0, 1, 0)
    xyj.add_atom("H3", h, -1, 0, 0)
    xyj.add_atom("H4", h, 0, -1, 0)
    xyj.add_atom("C", c, 0, 0, 0)
    mol = Xponge.Molecule(name="XYJ2")
    mol.add_residue(xyj)

    assert assign.atom_count == 3
    assert assign.bonds[0][1] == 1
    assert assign.bonds[1][2] == 1
    assert xyj.name == "XYJ"
    assert len(xyj.atoms) == 5
    assert mol.name == "XYJ2"
    assert mol.residue_count == 1


def test_legacy_assignrule_aliases_and_pure_string_path_match_old_usage():
    import Xponge
    from Xponge.assign import AssignRule

    rule = AssignRule("myrule", pure_string=True)

    @rule.Add_Rule("A", -1)
    def _a(i, assign):  # pylint: disable=unused-argument
        return True

    @rule.addRule("B", 1)
    def _b(i, assign):  # pylint: disable=unused-argument
        return True

    @rule.Set_Pre_Action
    def _pre(assign):
        assign.atoms[1] = "O"

    @rule.setPostAction
    def _post(assign):
        return assign

    assign = Xponge.Assign()
    assign.add_atom("H", 0, 0, 0)
    assign.add_atom("C", 1, 0, 0)
    assign.add_atom("O", 1, 0, 0)
    assign.add_bond(0, 1, 1)
    assign.add_bond(2, 1, 2)

    results = assign.determine_atom_type("myrule")
    assert results == ["B", "B", "B"]


def test_legacy_assign_delete_aliases_match_old_usage():
    import Xponge
    from io import StringIO

    s = StringIO(
        """12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
"""
    )
    ben = Xponge.get_assignment_from_xyz(s)
    ben.add_bond(0, 11, 1)
    ben.add_atom("O", 0, 0, 0)
    ben.add_bond(0, 12, 1)
    ben.Delete_Atom(12)

    assert ben.atom_count == 12
    ben.deleteBond(0, 11)
    assert ben.bond_count == 12
    assert 11 not in ben.bonds[0]
    assert 0 not in ben.bonds[11]


def test_legacy_assign_common_camelcase_aliases_exist():
    import Xponge

    assign = Xponge.Assign()
    names = [
        "Determine_Ring_And_Bond_Type",
        "Save_As_PDB",
        "Set_Charge",
        "Set_Charges",
        "Set_Coordinate",
        "Set_Formal_Charge",
        "Set_Atom_Type",
        "Add_Bond_Marker",
        "Has_Bond_Marker",
    ]
    for name in names:
        assert hasattr(assign, name), name


def test_assign_atom_type_accepts_atomtype_objects_and_property_setter():
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")

    assign.set_atom_type(0, hw)
    assign.Set_Atom_Type(1, hw)
    assign.atom_types = [hw, hw, "HW", "HW"]

    assert assign.atom_types == ["HW", "HW", "HW", "HW"]


def test_assign_to_residuetype_handles_duplicate_atom_names_like_legacy_xponge():
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")
    assign.atom_types = [hw, hw, "HW", "HW"]

    residue_type = assign.to_residuetype("TES")

    assert [atom.name for atom in residue_type.atoms] == ["O", "O1", "H", "H1"]


def test_save_sponge_input_accepts_residuetype_and_returns_molecule_like_legacy_xponge(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")
    assign.atom_types = [hw, hw, "HW", "HW"]
    residue_type = assign.to_residuetype("TES")

    monkeypatch.chdir(tmp_path)
    (tmp_path / "meta1d").mkdir()
    mol = Xponge.save_sponge_input(residue_type, "meta1d/test")

    assert mol is not None
    assert mol.residue_count == 1
    assert mol.atom_count == 4
    assert (tmp_path / "meta1d" / "test_coordinate.txt").is_file()


def test_legacy_metadynamics_workflow_front_half_runs_with_compat_layer(tmp_path, monkeypatch):
    from pathlib import Path

    import Xponge
    from Xponge.forcefield.base.angle_base import AngleType
    import Xponge.forcefield.amber.tip3p  # noqa: F401
    from Xponge.helper.cv import CVSystem

    monkeypatch.chdir(tmp_path)
    Path("meta1d").mkdir()

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")
    AngleType.New_From_String(
        """name       k   b
HW-HW-HW   50  1.7
"""
    )
    assign.atom_types = [hw, hw, hw, hw]
    tes = assign.to_residuetype("TES")
    mol = Xponge.save_sponge_input(tes, "meta1d/test")

    cv = CVSystem(mol)
    cv.add_cv_dihedral("torsion", mol.atoms[2], mol.atoms[0], mol.atoms[1], mol.atoms[3])
    cv.meta1d("torsion", CV_grid=10, CV_minimal=-3.142, CV_maximum=3.142, welltemp_factor=50, height=1, CV_sigma=0.5)
    cv.print("torsion")
    cv.output("meta1d/cv.txt")

    assert (tmp_path / "meta1d" / "test_coordinate.txt").is_file()
    assert (tmp_path / "meta1d" / "cv.txt").is_file()


def test_add_solvent_box_accepts_residuetype_like_legacy_fep_workflow():
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401

    assign = Xponge.get_assignment_from_smiles("C")
    residue_type = assign.to_residuetype("A")
    wat = Xponge.ResidueType.get_type("WAT")

    mol = Xponge.add_solvent_box(residue_type, wat, 8, n_solvent=3)

    assert mol.residue_count >= 2
    assert mol.residues[0].name == "A"


def test_legacy_fep_uncovalent_workflow_front_half_runs(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401
    import Xponge.forcefield.amber.gaff  # noqa: F401

    monkeypatch.chdir(tmp_path)

    assign_a = Xponge.get_assignment_from_smiles("C")
    assign_a.determine_atom_type("gaff")
    assign_a.calculate_charge("tpacm4")
    residue_a = assign_a.to_residuetype("A")
    wat = Xponge.ResidueType.get_type("WAT")
    solvated = Xponge.add_solvent_box(residue_a, wat, 10, n_solvent=3)
    Xponge.save_pdb(solvated, "test.pdb")
    Xponge.save_mol2(residue_a)

    assign_b = Xponge.get_assignment_from_smiles("CC")
    assign_b.determine_atom_type("gaff")
    assign_b.calculate_charge("tpacm4")
    Xponge.save_mol2(assign_b.to_residuetype("B"))

    assert (tmp_path / "test.pdb").is_file()
    assert (tmp_path / "A.mol2").is_file()
    assert (tmp_path / "B.mol2").is_file()


def test_legacy_fep_covalent_workflow_front_half_runs(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.tip3p  # noqa: F401
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    monkeypatch.chdir(tmp_path)

    ace = Xponge.ResidueType.get_type("ACE")
    wat = Xponge.ResidueType.get_type("WAT")
    ala = Xponge.ResidueType.get_type("ALA")
    gly = Xponge.ResidueType.get_type("GLY")
    nme = Xponge.ResidueType.get_type("NME")

    solvated = Xponge.add_solvent_box(ace + ala + nme, wat, 10, n_solvent=3)
    Xponge.save_pdb(solvated, "test.pdb")
    Xponge.save_mol2(ala, "A.mol2")
    Xponge.save_mol2(gly, "B.mol2")

    assert (tmp_path / "test.pdb").is_file()
    assert (tmp_path / "A.mol2").is_file()
    assert (tmp_path / "B.mol2").is_file()


def test_legacy_fep_helper_mutation_path_executes_and_exports_softcore(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    from Xponge.forcefield.special import fep

    monkeypatch.chdir(tmp_path)

    mol = Xponge.ACE + Xponge.ALA + Xponge.NME

    fep.Set_LJ_Type_B(mol, 0, "ZERO_LJ_ATOM")
    fep.Set_Subsys(mol, 0, 1)
    fep.Enable_LJ_Soft_Core(mol)
    fep.Save_Soft_Core_LJ(mol)
    outputs = Xponge.Save_SPONGE_Input(mol, "fep_soft")

    assert outputs is mol
    assert mol.atoms[0].lj_type_b == "ZERO_LJ_ATOM"
    assert mol.atoms[0].subsys == 1
    assert (tmp_path / "fep_soft_LJ_soft_core.txt").is_file()
    assert (tmp_path / "fep_soft_coordinate.txt").is_file()


def test_cvsystem_and_gb_expose_legacy_name_style_aliases():
    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401
    from Xponge.forcefield.special import gb
    from Xponge.helper.cv import CVSystem

    names = [
        "Add_Center",
        "AddCenter",
        "addCenter",
        "Add_CV_Position",
        "AddCVPosition",
        "addCvPosition",
        "Add_CV_Dihedral",
        "AddCVDihedral",
        "addCvDihedral",
        "Add_CV_RMSD",
        "AddCVRMSD",
        "addCvRMSD",
        "Print",
        "Steer",
        "Restrain",
        "Meta1D",
        "Meta_1D",
        "Output",
    ]
    for name in names:
        assert hasattr(CVSystem, name), name

    for name in ["set_gb_radius", "Set_GB_Radius", "SetGbRadius", "setGbRadius"]:
        assert hasattr(gb, name), name


def test_cvsystem_and_gb_camelcase_aliases_execute_real_workflow(tmp_path, monkeypatch):
    from pathlib import Path

    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401
    import Xponge.forcefield.amber.tip3p  # noqa: F401
    from Xponge.forcefield.special import gb
    from Xponge.helper.cv import CVSystem

    monkeypatch.chdir(tmp_path)
    mol = Xponge.NALA + Xponge.ALA + Xponge.CALA
    gb.SetGbRadius(mol)

    cv = CVSystem(mol)
    cv.Add_Center("c", "protein")
    cv.Add_CV_Position("x", "c", "x", False)
    cv.Print("x")
    cv.Steer("x", 2)
    Path("cv").mkdir()
    cv.Output("cv/cv.txt", folder="cv")

    text = (tmp_path / "cv" / "cv.txt").read_text()
    assert "vatom_type = center" in text
    assert "CV_type = position_x" in text
    assert "steer" in text


def test_cvsystem_restrain_and_meta1d_camelcase_aliases_execute(tmp_path, monkeypatch):
    from pathlib import Path

    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401
    from Xponge.helper.cv import CVSystem

    monkeypatch.chdir(tmp_path)
    mol = Xponge.NALA + Xponge.ALA + Xponge.CALA

    cv = CVSystem(mol)
    cv.Add_CV_Dihedral("torsion", mol.atoms[2], mol.atoms[0], mol.atoms[1], mol.atoms[3])
    cv.Restrain("torsion", 200, "ref")
    cv.Meta1D(
        "torsion",
        CV_grid=10,
        CV_minimal=-3.14,
        CV_maximum=3.14,
        CV_sigma=0.5,
        welltemp_factor=20,
        height=1,
    )
    Path("cv").mkdir()
    cv.Output("cv/cv.txt", folder="cv")

    text = (tmp_path / "cv" / "cv.txt").read_text()
    assert "restrain" in text
    assert "reference = ref" in text
    assert "metad" in text
    assert "CV_period" in text


def test_cvsystem_rmsd_camelcase_alias_executes(tmp_path, monkeypatch):
    from pathlib import Path

    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
    from Xponge.helper.cv import CVSystem

    monkeypatch.chdir(tmp_path)
    mol = Xponge.ACE + Xponge.ALA * 2 + Xponge.NME

    cv = CVSystem(mol)
    cv.Add_CV_RMSD("RMSD", "protein and backbone")
    cv.Print("RMSD")
    Path("cv").mkdir()
    cv.Output("cv/cv.txt", folder="cv")

    assert (tmp_path / "cv" / "RMSD_atom.txt").is_file()
    assert (tmp_path / "cv" / "RMSD_coordinate.txt").is_file()
    text = (tmp_path / "cv" / "cv.txt").read_text()
    assert "CV_type = rmsd" in text
    assert "CV = RMSD" in text


def test_top_level_load_alias_variants_cover_legacy_name_style_examples():
    import Xponge

    names = [
        "Load_PDB",
        "Load_Pdb",
        "LoadPDB",
        "LoadPdb",
        "Load_GRO",
        "Load_Gro",
        "LoadGRO",
        "LoadGro",
        "Load_Mol2",
        "LoadMol2",
    ]
    for name in names:
        assert hasattr(Xponge, name), name


def test_xponge_package_alias_load_pdb_variants_support_real_io_options():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    pdb_text = r"""
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
    from io import StringIO

    p1 = Xponge.Load_PDB(StringIO(pdb_text), ignore_hydrogen=True, unterminal_residues=[2, 3])
    p2 = Xponge.LoadPdb(StringIO(pdb_text), ignore_hydrogen=True, unterminal_residues=["A:2"])

    assert [res.name for res in p1.residues] == ["VAL", "TRP"]
    assert [res.name for res in p2.residues] == ["VAL", "CTRP"]


def test_xponge_package_alias_save_pdb_variant_handles_write_cryst1(tmp_path):
    import Xponge
    import Xponge.forcefield.base.mass_base  # noqa: F401
    import Xponge.forcefield.base.charge_base  # noqa: F401

    Xponge.AtomType.New_From_String(
        """
name  mass   charge[e]
XLP   12.00  0.00
"""
    )
    residue_type = Xponge.ResidueType(name="GLN__HE21_H__HE22_H")
    residue_type.add_atom("CA", "XLP", 72.610, 47.770, 57.850)
    mol = Xponge.Molecule(name="LONGRES")
    mol.add_residue(residue_type)
    mol.box_length = [91.0, 91.0, 91.0]

    outfile = tmp_path / "longres_alias.pdb"
    Xponge.SavePDB(mol, str(outfile), write_cryst1=False)

    lines = outfile.read_text().splitlines()
    atom_line = next(line for line in lines if line.startswith("ATOM"))
    assert not lines[0].startswith("CRYST1")
    assert atom_line[17:20] == "GLN"
    assert float(atom_line[30:38]) == 72.610


def test_xponge_package_alias_loadmol2_and_savemol2_support_old_io_shape(tmp_path, monkeypatch):
    import Xponge
    from io import StringIO

    mol2_contents = r"""
@<TRIPOS>MOLECULE
ASN
 12 12 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    C  -1.2131  -0.6884   0.0000   ca            1      ASN   0.000000
     2   C1  -1.2028   0.7064   0.0001   ca            1      ASN   0.000000
     3   C2  -0.0103  -1.3948   0.0000   ca            1      ASN   0.000000
     4   C3   0.0104   1.3948  -0.0001   ca            1      ASN   0.000000
     5   C4   1.2028  -0.7063   0.0000   ca            1      ASN   0.000000
     6   C5   1.2131   0.6884   0.0000   ca            1      ASN   0.000000
     7    H  -2.1577  -1.2244   0.0000   ha            1      ASN   0.000000
     8   H1  -2.1393   1.2564   0.0001   ha            1      ASN   0.000000
     9   H2  -0.0184  -2.4809  -0.0001   ha            1      ASN   0.000000
    10   H3   0.0184   2.4808   0.0000   ha            1      ASN   0.000000
    11   H4   2.1394  -1.2563   0.0001   ha            1      ASN   0.000000
    12   H5   2.1577   1.2245   0.0000   ha            1      ASN   0.000000
@<TRIPOS>BOND
     1      1      2 ar
     2      1      3 ar
     3      1      7 1
     4      2      4 ar
     5      2      8 1
     6      3      5 ar
     7      3      9 1
     8      4      6 ar
     9      4     10 1
    10      5      6 ar
    11      5     11 1
    12      6     12 1
@<TRIPOS>SUBSTRUCTURE
    1      ASN      1 ****               0 ****  ****
"""

    ben = Xponge.LoadMol2(StringIO(mol2_contents.replace("ASN", "BBB")), ignore_atom_type=True)
    assert Xponge.ResidueType.get_type("BBB").name == "BBB"

    import Xponge.forcefield.amber.gaff  # noqa: F401

    Xponge.Load_Mol2(StringIO(mol2_contents.replace("ASN", "BEN")), as_template=True)
    assert Xponge.ResidueType.get_type("BEN").name == "BEN"

    outfile = tmp_path / "ben_alias.mol2"
    Xponge.SaveMol2(ben, str(outfile))
    assert outfile.is_file()
    assert "@<TRIPOS>MOLECULE" in outfile.read_text()

    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    monkeypatch.chdir(tmp_path)
    Xponge.save_mol2(Xponge.ResidueType.get_type("ALA"))
    assert (tmp_path / "ALA.mol2").is_file()


def test_xponge_package_alias_save_gro_supports_legacy_io_shape(tmp_path):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.NALA + Xponge.ALA + Xponge.CALA
    outfile = tmp_path / "legacy.gro"
    Xponge.save_gro(mol, str(outfile))

    assert outfile.is_file()
    lines = outfile.read_text().splitlines()
    assert lines[0] == "Generated By Xponge"
    assert int(lines[1].strip()) == mol.atom_count


def test_xponge_package_alias_save_gro_camelcase_and_instance_shapes_work(tmp_path):
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.NALA + Xponge.ALA + Xponge.CALA
    camel = tmp_path / "camel.gro"
    inst = tmp_path / "instance.gro"

    Xponge.Save_GRO(mol, str(camel))
    mol.Save_GRO(str(inst))

    assert camel.is_file()
    assert inst.is_file()
    assert camel.read_text().splitlines()[0] == "Generated By Xponge"
    assert inst.read_text().splitlines()[0] == "Generated By Xponge"


def test_forcefield_import_syncs_template_globals_into_xponge_modules(tmp_path, monkeypatch):
    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401
    from Xponge.forcefield.special import gb

    assert hasattr(Xponge, "NALA")
    assert hasattr(Xponge, "ALA")
    assert hasattr(Xponge, "CALA")

    monkeypatch.chdir(tmp_path)
    mol = Xponge.NALA + Xponge.ALA * 2 + Xponge.CALA
    gb.set_gb_radius(mol)
    Xponge.Save_PDB(mol, "ALA.pdb")
    Xponge.Save_SPONGE_Input(mol, "ALA")

    assert (tmp_path / "ALA.pdb").is_file()
    assert (tmp_path / "ALA_gb.txt").is_file()


def test_xponge_import_injects_legacy_top_level_names_into_main_namespace(tmp_path, monkeypatch):
    import __main__

    import Xponge
    import Xponge.forcefield.amber.ff19sb  # noqa: F401
    from Xponge.forcefield.special import gb

    assert "Save_PDB" in __main__.__dict__
    assert "Save_SPONGE_Input" in __main__.__dict__
    assert "NALA" in __main__.__dict__
    assert "ALA" in __main__.__dict__
    assert "CALA" in __main__.__dict__

    monkeypatch.chdir(tmp_path)
    __main__.__dict__["gb"] = gb
    exec(
        "mol = NALA + ALA * 2 + CALA\n"
        "gb.Set_GB_Radius(mol)\n"
        "Save_PDB(mol, 'ALA.pdb')\n"
        "Save_SPONGE_Input(mol, 'ALA')\n",
        __main__.__dict__,
    )

    assert (tmp_path / "ALA.pdb").is_file()
    assert (tmp_path / "ALA_gb.txt").is_file()


def test_helper_base_class_shapes_are_minimally_compatible():
    import Xponge
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.get_template_molecule("ALA")
    assert isinstance(mol, Xponge.AbstractMolecule)
    link = Xponge.ResidueLink(mol.atoms[0], mol.atoms[1])
    assert isinstance(link, Xponge.Entity)


def test_xpongecpp_legacy_alias_names_remain_callable():
    import XpongeCPP as Xponge

    assert Xponge.LoadPDB is Xponge.load_pdb
    assert Xponge.SavePDB is Xponge.Save_PDB
    assert Xponge.GetAssignmentFromPDB is Xponge.get_assignment_from_pdb
    assert Xponge.AddMolecule is Xponge.Add_Molecule
    assert Xponge.AddSolventBox is Xponge.Add_Solvent_Box
    assert Xponge.SetBoxPadding is Xponge.Set_Box_Padding


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


def test_helper_compat_exposes_globalsetting_runtime_shape():
    import Xponge
    from Xponge.helper import Xopen, remove_real_global_variable, set_real_global_variable, source

    assert Xponge.GlobalSetting.purpose == "academic"
    assert hasattr(Xponge.GlobalSetting, "logger")
    assert hasattr(Xponge.GlobalSetting.logger, "handlers")
    assert callable(Xopen)
    assert callable(source)
    namespace = {}
    set_real_global_variable("TMP", 1, namespace=namespace)
    remove_real_global_variable("TMP", namespace=namespace)
    assert "TMP" not in namespace


def test_globalsetting_legacy_method_aliases_exist():
    import Xponge

    Xponge.GlobalSetting.Add_GMX_Include_Path("/tmp/demo")
    assert "/tmp/demo" in Xponge.GlobalSetting.GMXIncludePaths
    Xponge.GlobalSetting.add_pdb_residue_alias_mapping("H2O", "WAT")
    Xponge.GlobalSetting.set_invisible_bonded_forces(["bond"])

    decorator = Xponge.GlobalSetting.Add_Unit_Transfer_Function(object)
    wrapped = decorator(lambda self: self)

    assert callable(wrapped)


def test_helper_source_supports_forcefield_import():
    from Xponge.helper import source

    module = source("Xponge.forcefield.amber.ff14sb", into_global=False)

    assert module is not None
    assert module.__name__.endswith("Xponge.forcefield.amber.ff14sb")


def test_helper_exports_bonded_force_generators_and_data_amber_imports():
    import Xponge
    import XpongeCPP
    import XpongeCPP.data.amber as amber_data
    import Xponge.forcefield.amber.rsff2c as rsff2c
    from Xponge.helper import Generate_New_Bonded_Force_Type
    from XpongeCPP.helper import generate_new_bonded_force_type

    assert amber_data is not None
    assert rsff2c is not None
    assert callable(Generate_New_Bonded_Force_Type)
    assert callable(generate_new_bonded_force_type)
    assert callable(Xponge.Generate_New_Bonded_Force_Type)
    assert callable(XpongeCPP.generate_new_bonded_force_type)


def test_public_compat_module_exports_runtime_installers():
    import XpongeCPP.compat as compat

    assert callable(compat.install_legacy_assign_patches)
    assert callable(compat.Install_Legacy_Assign_Patches)
    assert callable(compat.install_legacy_runtime_patches)
    assert callable(compat.Install_Legacy_Runtime_Patches)
    assert callable(compat.get_legacy_residue_links_override)
    assert callable(compat.Get_Legacy_Residue_Links_Override)


def test_save_pdb_preserves_python_side_residue_link_overrides(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.ResidueType.get_type("ALA") * 2
    atom_a = mol.residues[0].atoms[0]
    atom_b = mol.residues[1].atoms[0]
    mol.residue_links = [(atom_a, atom_b)]

    target = tmp_path / "links.pdb"
    Xponge.Save_PDB(mol, target)
    text = target.read_text()

    assert "CONECT" in text


def test_forcefield_base_first_wave_minimal_ffbase_compatibility():
    import numpy as np
    import Xponge
    import Xponge.forcefield.base.mass_base  # noqa: F401
    import Xponge.forcefield.base.charge_base  # noqa: F401
    from Xponge.forcefield.base.bond_base import BondType, _gmx_parser
    from Xponge.forcefield.base.cmap_base import CMapType
    from Xponge.forcefield.base.dihedral_base import ImproperType, ProperType
    from Xponge.forcefield.base.lj_base import LJType
    from Xponge.forcefield.base.nb14_base import NB14Type
    from Xponge.forcefield.base.virtual_atom_base import VirtualType2

    LJType.New_From_String(
        r"""
name      A         B
AG-AG     2021000   6072
AL-AL     1577000   5035
AU-AU     2307000   6987
"""
    )
    LJType.New_From_String(
        """
name    epsilon   rmin
ag-ag     4.56    1.4775
al-al     4.02    1.4625
au-au     5.29    1.4755
"""
    )
    for name in ["Ag-Ag", "Al-Al", "Au-Au"]:
        er = LJType.get_type(name.lower())
        ab = LJType.get_type(name.upper())
        assert abs(er.epsilon - ab.epsilon) < 0.01
        assert abs(er.rmin - ab.rmin) < 0.01

    LJType.New_From_String(
        r"""
name    epsilon[eV]   sigma[nm]
y-y     0.0017345     0.32
"""
    )
    yy = LJType.get_type("y-y")
    assert abs(yy.epsilon - 0.03999851) < 0.01
    assert abs(yy.rmin - np.power(2, 1 / 6) * 3.2 / 2) < 0.01

    Xponge.AtomType.New_From_String(
        r"""
name    mass    charge[e]   LJtype
H       1.008   1.000       HW
"""
    )
    h = Xponge.AtomType.get_type("H")
    assert h.mass == 1.008
    assert h.__dict__["charge[e]"] == 1.0
    assert h.LJtype == "HW"

    BondType.New_From_String(
        """
name k b
XTBOND_A-XTBOND_B 100.0 1.5
XTBOND_D-XTBOND_C 200.0 1.6
"""
    )

    class DummyAtomType:
        def __init__(self, name):
            self.name = name

    class DummyAtom:
        def __init__(self, type_):
            self.type = type_

    class DummyMol:
        def __init__(self):
            self.forces = []

        def add_bonded_force(self, force):
            self.forces.append(force)

    mol1 = DummyMol()
    stat1 = {1: DummyAtom(DummyAtomType("XTBOND_A")), 2: DummyAtom(DummyAtomType("XTBOND_B"))}
    _gmx_parser(["1", "2", "1"], mol1, stat1)
    assert len(mol1.forces) == 1

    mol2 = DummyMol()
    stat2 = {1: DummyAtom(DummyAtomType("XTBOND_C")), 2: DummyAtom(DummyAtomType("XTBOND_D"))}
    _gmx_parser(["1", "2", "1"], mol2, stat2)
    assert len(mol2.forces) == 1

    mol3 = DummyMol()
    stat3 = {1: DummyAtom(DummyAtomType("XTBOND_MISSING")), 2: DummyAtom(DummyAtomType("XTBOND_ALSO_MISSING"))}
    try:
        _gmx_parser(["1", "2", "1"], mol3, stat3)
    except KeyError:
        pass
    else:  # pragma: no cover
        raise AssertionError("Expected KeyError for missing bond type")

    ProperType.Clear_Type()
    ImproperType.Clear_Type()
    ProperType.New_From_String(
        """
name k phi0 periodicity
CT-CT-CT-CT 0.156 3.1415926 3
"""
    )
    ImproperType.New_From_String(
        """
name k phi0 periodicity
CA-CA-C-O 1.1 3.1415926 2
"""
    )
    proper = ProperType.get_type("ct-ct-ct-ct")
    improper = ImproperType.get_type("CA-CA-C-O")
    assert proper.k == pytest.approx(0.156)
    assert proper.periodicity == 3
    assert improper.k == pytest.approx(1.1)
    assert improper.periodicity == 2
    same_force = ImproperType.same_force("A-B-C-D")
    assert "A-B-C-D" in same_force
    assert "B-A-C-D" in same_force

    NB14Type.Clear_Type()
    NB14Type.New_From_String(
        """
name kLJ kee
X-X 0.5 0.833333
CT-H1 0.25 0.5
"""
    )
    xx = NB14Type.get_type("x-x")
    cth1 = NB14Type.Get_Type("CT-H1")
    assert xx.kLJ == pytest.approx(0.5)
    assert xx.kee == pytest.approx(0.833333)
    assert cth1.kLJ == pytest.approx(0.25)
    assert cth1.kee == pytest.approx(0.5)

    VirtualType2.Clear_Type()
    VirtualType2.New_From_String(
        """
name atom0 atom1 atom2 k1 k2
EP -3 -2 -1 0.1066413 0.1066413
"""
    )
    ep = VirtualType2.get_type("EP")
    assert ep.atom0 == -3
    assert ep.atom1 == -2
    assert ep.atom2 == -1
    assert ep.k1 == pytest.approx(0.1066413)
    assert ep.k2 == pytest.approx(0.1066413)
    assert Xponge.GlobalSetting.VirtualAtomTypes["vatom2"] == 3

    CMapType.Clear_Type()
    CMapType.New_From_Dict(
        {
            "ALA@C-N-CA-C-N": {
                "resolution": 2,
                "parameters": [1.0, 2.0, 3.0, 4.0],
            }
        }
    )
    cmap = CMapType.get_type("ALA@C-N-CA-C-N")
    assert cmap.resolution == 2
    assert cmap.parameters == [1.0, 2.0, 3.0, 4.0]
    assert CMapType.same_force("A-B-C-D-E") == ["A-B-C-D-E"]
