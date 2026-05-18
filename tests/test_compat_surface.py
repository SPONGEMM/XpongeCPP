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


def test_public_compat_bootstrap_entrypoint_is_available_without_hidden_auto_api():
    import XpongeCPP.compat as compat

    assert callable(compat.install_legacy_bootstrap)
    assert callable(compat.enable_legacy_namespace)


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
    from XpongeCPP.mdrun import _BIN_PATH_FILE

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
    assert "src/XpongeCPP/BIN_PATH.dat" not in str(_BIN_PATH_FILE)


def test_old_xponge_init_export_surface_has_first_wave_compatibility_shape():
    import Xponge
    from Xponge import build_bonded_force, get_mindsponge_system_energy, set_global_alternative_names

    assert callable(set_global_alternative_names)
    assert callable(build_bonded_force)
    assert callable(get_mindsponge_system_energy)
    assert hasattr(Xponge, "Atom")
    assert hasattr(Xponge, "Residue")
