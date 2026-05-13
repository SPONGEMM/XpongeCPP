import math
from io import StringIO

import numpy as np
import pytest

import XpongeCPP as Xponge


def _load_chain():
    return Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
CHAIN
4 3 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
2 A2 1.0 0.0 0.0 C 1 LIG 0.0
3 A3 2.0 1.0 0.0 C 1 LIG 0.0
4 A4 3.0 1.0 1.0 C 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
"""
        )
    )


def _distance(atom1, atom2):
    return math.sqrt((atom1.x - atom2.x) ** 2 + (atom1.y - atom2.y) ** 2 + (atom1.z - atom2.z) ** 2)


def _angle(atom1, atom2, atom3):
    v1 = np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z])
    v2 = np.array([atom3.x - atom2.x, atom3.y - atom2.y, atom3.z - atom2.z])
    cosine = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
    cosine = max(-1.0, min(1.0, float(cosine)))
    return math.acos(cosine)


def _dihedral(atom1, atom2, atom3, atom4):
    p0 = np.array([atom1.x, atom1.y, atom1.z], dtype=float)
    p1 = np.array([atom2.x, atom2.y, atom2.z], dtype=float)
    p2 = np.array([atom3.x, atom3.y, atom3.z], dtype=float)
    p3 = np.array([atom4.x, atom4.y, atom4.z], dtype=float)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.atan2(y, x)


def test_impose_bond_updates_target_distance():
    mol = _load_chain()
    atom1, atom2 = mol.residues[0].atoms[:2]

    Xponge.impose_bond(mol, atom1, atom2, 2.0)

    assert math.isclose(_distance(atom1, atom2), 2.0, rel_tol=0.0, abs_tol=1e-6)


def test_impose_angle_updates_target_angle():
    mol = _load_chain()
    atom1, atom2, atom3 = mol.residues[0].atoms[:3]

    Xponge.impose_angle(mol, atom1, atom2, atom3, math.pi / 2)

    assert math.isclose(_angle(atom1, atom2, atom3), math.pi / 2, rel_tol=0.0, abs_tol=1e-6)


def test_impose_dihedral_updates_target_dihedral():
    mol = _load_chain()
    atom1, atom2, atom3, atom4 = mol.residues[0].atoms

    Xponge.impose_dihedral(mol, atom1, atom2, atom3, atom4, math.pi / 3)

    assert math.isclose(_dihedral(atom1, atom2, atom3, atom4), math.pi / 3, rel_tol=0.0, abs_tol=1e-6)


def test_main_axis_rotate_aligns_long_axis_to_requested_direction():
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
LINE
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
2 A2 5.0 0.0 0.0 C 1 LIG 0.0
3 A3 10.0 0.0 0.0 C 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""
        )
    )

    Xponge.main_axis_rotate(mol, direction_long=[0, 0, 1])

    xs = [atom.x for atom in mol.residues[0].atoms]
    zs = [atom.z for atom in mol.residues[0].atoms]
    assert max(xs) - min(xs) < 1e-6
    assert max(zs) - min(zs) > 9.9


def test_get_peptide_from_sequence_uses_terminal_templates():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.get_peptide_from_sequence("AG")

    assert [res.name for res in mol.residues] == ["NALA", "CGLY"]


def test_h_mass_repartition_moves_mass_from_heavy_atom_to_hydrogen():
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
CH
2 1 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C 1 LIG 0.0
2 H1 1.0 0.0 0.0 H 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
"""
        )
    )

    atom_c, atom_h = mol.residues[0].atoms
    Xponge.h_mass_repartition(mol, repartition_mass=1.1, repartition_rate=3)

    assert math.isclose(atom_h.mass, 3.024, rel_tol=0.0, abs_tol=1e-6)
    assert math.isclose(atom_c.mass, 9.994, rel_tol=0.0, abs_tol=1e-6)


def test_solvent_replace_replaces_selected_residues_and_sorts():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    solute = Xponge.get_peptide_from_sequence("AG")
    water = Xponge.get_template_molecule("WAT")
    Xponge.Add_Solvent_Box(solute, water, 4.0, tolerance=2.5, n_solvent=6, seed=20260509)

    before = solute.residue_count
    Xponge.solvent_replace(solute, water, {Xponge.get_template_molecule("NA"): 1, Xponge.get_template_molecule("CL"): 1})

    counts = solute.residue_counts()
    assert counts["NA"] == 1
    assert counts["CL"] == 1
    assert counts["WAT"] == 4
    assert solute.residue_count == before
    assert [res.name for res in solute.residues[:2]] == ["NALA", "CGLY"]
    assert solute.validate()


def test_solvent_replace_supports_mixed_replacements_with_deterministic_seed():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    def build_system():
        mol = Xponge.get_peptide_from_sequence("AG")
        Xponge.Add_Solvent_Box(mol, Xponge.get_template_molecule("WAT"), 4.0, tolerance=2.5, n_solvent=8, seed=20260509)
        return mol

    np.random.seed(20260510)
    mol_a = build_system()
    Xponge.solvent_replace(
        mol_a,
        Xponge.get_template_molecule("WAT"),
        {
            Xponge.get_template_molecule("NA"): 2,
            Xponge.get_template_molecule("CL"): 1,
            Xponge.get_template_molecule("K"): 1,
        },
    )

    np.random.seed(20260510)
    mol_b = build_system()
    Xponge.solvent_replace(
        mol_b,
        Xponge.get_template_molecule("WAT"),
        {
            Xponge.get_template_molecule("NA"): 2,
            Xponge.get_template_molecule("CL"): 1,
            Xponge.get_template_molecule("K"): 1,
        },
    )

    assert [res.name for res in mol_a.residues] == [res.name for res in mol_b.residues]
    assert mol_a.residue_counts() == {"NALA": 1, "CGLY": 1, "NA": 2, "CL": 1, "K": 1, "WAT": 4}
    assert [res.name for res in mol_a.residues[:2]] == ["NALA", "CGLY"]
    assert mol_a.validate()


def test_sort_atoms_by_reorders_residue_atoms_to_match_template():
    template = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
TEMPLATE
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
2 A2 1.0 0.0 0.0 O 1 LIG 0.0
3 A3 2.0 0.0 0.0 H 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""
        )
    )
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
MOL
3 2 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A3 2.0 0.0 0.0 H 1 LIG 0.0
2 A1 0.0 0.0 0.0 C 1 LIG 0.0
3 A2 1.0 0.0 0.0 O 1 LIG 0.0
@<TRIPOS>BOND
1 2 3 1
2 3 1 1
"""
        )
    )

    Xponge.sort_atoms_by(mol, template)

    assert [atom.name for atom in mol.residues[0].atoms] == ["A1", "A2", "A3"]
    assert sorted(tuple(sorted(bond)) for bond in mol.explicit_bonds) == [(0, 1), (1, 2)]
    assert mol.validate()


def test_sort_atoms_by_reorders_multiple_residues_independently():
    template = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
TEMPLATE2
6 4 2
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
2 A2 1.0 0.0 0.0 O 1 LIG 0.0
3 A3 2.0 0.0 0.0 H 1 LIG 0.0
4 B1 4.0 0.0 0.0 C 2 LIG 0.0
5 B2 5.0 0.0 0.0 O 2 LIG 0.0
6 B3 6.0 0.0 0.0 H 2 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 4 5 1
4 5 6 1
"""
        )
    )
    mol = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
MOL2
6 4 2
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A3 2.0 0.0 0.0 H 1 LIG 0.0
2 A1 0.0 0.0 0.0 C 1 LIG 0.0
3 A2 1.0 0.0 0.0 O 1 LIG 0.0
4 B2 5.0 0.0 0.0 O 2 LIG 0.0
5 B3 6.0 0.0 0.0 H 2 LIG 0.0
6 B1 4.0 0.0 0.0 C 2 LIG 0.0
@<TRIPOS>BOND
1 2 3 1
2 3 1 1
3 6 4 1
4 4 5 1
"""
        )
    )

    Xponge.sort_atoms_by(mol, template)

    assert [atom.name for atom in mol.residues[0].atoms] == ["A1", "A2", "A3"]
    assert [atom.name for atom in mol.residues[1].atoms] == ["B1", "B2", "B3"]
    assert mol.validate()


def test_original_process_workflow_shape_is_reproducible_without_external_sponge(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.get_peptide_from_sequence("AAAAAAAAAA")
    Xponge.main_axis_rotate(mol)
    Xponge.add_solvent_box(mol, Xponge.get_template_molecule("WAT"), 20.0, seed=20260509)
    Xponge.solvent_replace(
        mol,
        Xponge.get_template_molecule("WAT"),
        {
            Xponge.get_template_molecule("NA"): 5,
            Xponge.get_template_molecule("CL"): 5,
        },
    )
    Xponge.h_mass_repartition(mol)

    outputs = Xponge.save_sponge_input(mol, prefix="ala", dirname=str(tmp_path))

    counts = mol.residue_counts()
    assert counts
    assert counts["NA"] == 5
    assert counts["CL"] == 5
    assert counts["WAT"] > 0
    assert {"bond", "angle", "dihedral", "exclude", "nb14"}.issubset(outputs)
    assert mol.validate()


def test_optimize_reports_missing_engine_explicitly(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _name: None)

    with pytest.raises(RuntimeError, match="SPONGE executable is required"):
        Xponge.optimize(_load_chain())


def test_optimize_respects_arguments_and_reports_engine_failure(monkeypatch, tmp_path):
    mol = _load_chain()
    calls = {}
    original_box_length = list(mol.box_length)
    original_box_angle = list(mol.box_angle)

    def fake_save(molecule, prefix="", dirname=""):
        calls["box_length_during_save"] = list(molecule.box_length)
        calls["box_angle_during_save"] = list(molecule.box_angle)
        return {"bond"}

    def fake_run(command, check=False, text=True, capture_output=True):
        calls["command"] = command
        class Result:
            returncode = 1
            stdout = ""
            stderr = "engine failed"
        return Result()

    monkeypatch.setattr("shutil.which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(Xponge, "Save_SPONGE_Input", fake_save)
    monkeypatch.setattr("subprocess.run", fake_run)

    with pytest.raises(RuntimeError, match="engine failed"):
        Xponge.optimize(
            mol,
            step=123,
            dt=2e-8,
            pbc=False,
            only_bad_coordinate=False,
            extra_commands={"custom_flag": 5},
        )

    assert calls["box_length_during_save"] == [999.0, 999.0, 999.0]
    assert calls["box_angle_during_save"] == [90.0, 90.0, 90.0]
    assert mol.box_length == original_box_length
    assert mol.box_angle == original_box_angle
    assert calls["command"][0] == "SPONGE_NOPBC"
    assert calls["command"][1] == "-mdin"
    assert isinstance(calls["command"][2], str)
    assert calls["command"][3:] == ["-dt", "2e-08"]


def test_optimize_runs_real_engine_when_available(monkeypatch):
    import os
    import shutil

    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    executable = shutil.which("SPONGE_NOPBC")
    pbc = False
    if executable is None:
        executable = shutil.which("SPONGE")
        pbc = True
    if executable is None:
        pytest.skip("local SPONGE executables are not available")

    mol = Xponge.get_peptide_from_sequence("AG")
    if pbc:
        mol.box_length = [20.0, 20.0, 20.0]
        mol.box_angle = [90.0, 90.0, 90.0]

    original_load_coordinate = Xponge.load_coordinate
    loaded = {}

    def tracking_load_coordinate(path, molecule):
        loaded["path"] = path
        return original_load_coordinate(path, molecule)

    monkeypatch.setattr(Xponge, "load_coordinate", tracking_load_coordinate)
    Xponge.optimize(mol, step=1, only_bad_coordinate=False, pbc=pbc)

    assert loaded["path"].endswith("_coordinate.txt")
    assert os.path.exists(loaded["path"])
    assert all(math.isfinite(atom.x + atom.y + atom.z) for atom in mol.atoms)
    assert mol.validate()


def test_lattice_create_populates_block_region_with_simple_cubic_basis():
    basis = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
ATOM
1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
"""
        )
    )

    box = Xponge.BlockRegion(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, boundary=True)
    lattice = Xponge.Lattice("sc", basis_molecule=basis, scale=1.0)

    mol = lattice.create(box, box)

    assert mol.atom_count == 8
    assert mol.residue_count == 8
    assert list(mol.box_length) == [2.0, 2.0, 2.0]


def test_lattice_create_supports_hcp_prism_region():
    basis = Xponge.load_mol2(
        StringIO(
            """@<TRIPOS>MOLECULE
ATOM
1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 C 1 LIG 0.0
"""
        )
    )

    prism = Xponge.PrismRegion(
        0.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        1.0, math.sqrt(3.0), 0.0,
        0.0, 0.0, 2.0,
        boundary=True,
    )
    lattice = Xponge.Lattice("hcp", basis_molecule=basis, scale=1.0)

    mol = lattice.create(prism, prism)

    assert mol.atom_count > 0
    assert mol.residue_count == mol.atom_count
    assert mol.box_length[0] > 0
    assert mol.box_angle[2] == pytest.approx(60.0)
    assert mol.validate()
