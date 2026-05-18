import XpongeCPP as Xponge
import pytest


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


def test_xponge_main_namespace_injection_can_be_disabled_by_policy(monkeypatch):
    import __main__
    import importlib
    import sys

    monkeypatch.setenv("XPONGECPP_LEGACY_MAIN_NAMESPACE", "0")
    for name in ("Save_PDB", "Save_SPONGE_Input", "NALA", "ALA", "CALA"):
        __main__.__dict__.pop(name, None)

    sys.modules.pop("Xponge", None)
    import Xponge as Xponge_alias  # noqa: F401

    importlib.reload(sys.modules["Xponge"])

    assert "Save_PDB" not in __main__.__dict__
    assert "Save_SPONGE_Input" not in __main__.__dict__
    assert "NALA" not in __main__.__dict__


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


def test_internal_compat_package_no_longer_reexports_public_helper_surface():
    import XpongeCPP._compat as internal_compat
    import XpongeCPP.compat as public_compat

    assert "api" in internal_compat.__all__
    assert "runtime" in internal_compat.__all__
    assert not hasattr(internal_compat, "enable_legacy_namespace")
    assert not hasattr(internal_compat, "install_legacy_runtime_patches")
    assert callable(public_compat.enable_legacy_namespace)
    assert callable(public_compat.install_legacy_runtime_patches)


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
