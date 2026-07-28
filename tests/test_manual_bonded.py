from dataclasses import replace
from io import StringIO
import os
from pathlib import Path
import subprocess

import pytest
import XpongeCPP as Xponge
from XpongeCPP.io_bundle import convert_bundle_to_legacy
from XpongeCPP.metal_assignment import (
    AtomParameterUpdate,
    ElectronicState,
    LJParameter,
    MANUAL_BONDED_SOURCE,
    ManualAngle,
    ManualBond,
    MetalAssignmentValidationError,
    MetalSite,
    apply_manual_bonded_assignment,
    molecule_input_hash,
    prepare_manual_bonded_assignment,
)


METAL_ANGLE_PDB = """\
HETATM    1 ZN   MZN A   1       0.000   0.000   0.000  1.00  0.00          Zn
HETATM    2 N1   LIG A   2       2.100   0.000   0.000  1.00  0.00           N
HETATM    3 O1   OLG A   3       0.000   2.000   0.000  1.00  0.00           O
END
"""

SAME_RESIDUE_PDB = """\
HETATM    1 ZN   CMP A   1       0.000   0.000   0.000  1.00  0.00          Zn
HETATM    2 N1   CMP A   1       2.100   0.000   0.000  1.00  0.00           N
HETATM    3 O1   CMP A   1       0.000   2.000   0.000  1.00  0.00           O
END
"""


def _load(text=METAL_ANGLE_PDB):
    return Xponge.load_pdb(StringIO(text))


def _prepare(molecule, *, angles=None):
    if angles is None:
        angles = (ManualAngle((1, 0, 2), 45.0, 1.5707963267948966),)
    return prepare_manual_bonded_assignment(
        molecule,
        metal_sites=(MetalSite(0, "Zn", 2),),
        electronic_state=ElectronicState(0, 1),
        coordination_edges=((0, 1), (0, 2)),
        bonds=(
            ManualBond((0, 1), 120.0, 2.1),
            ManualBond((0, 2), 140.0, 2.0),
        ),
        angles=angles,
    )


def test_manual_bonded_plan_is_hash_closed_and_side_effect_free():
    molecule = _load()
    input_hash = molecule_input_hash(molecule)

    plan = _prepare(molecule)

    assert plan.plan_hash == plan.computed_hash()
    assert plan.reference.artifact_hash == plan.reference.computed_hash()
    assert plan.assignment_plan.parameter_overlay.parameter_source == (
        MANUAL_BONDED_SOURCE
    )
    assert molecule_input_hash(molecule) == input_hash
    assert molecule.residue_links == []


def test_manual_bonded_applies_bonds_angles_and_preserves_parent():
    molecule = _load()

    result = apply_manual_bonded_assignment(molecule, _prepare(molecule))

    assert result.molecule is not molecule
    assert result.molecule.residue_links == [[0, 1], [0, 2]]
    assert molecule.residue_links == []


def test_manual_bonded_requires_complete_bond_and_angle_coverage():
    molecule = _load()

    with pytest.raises(MetalAssignmentValidationError) as exc:
        prepare_manual_bonded_assignment(
            molecule,
            metal_sites=(MetalSite(0, "Zn", 2),),
            electronic_state=ElectronicState(0, 1),
            coordination_edges=((0, 1), (0, 2)),
            bonds=(ManualBond((0, 1), 120.0, 2.1),),
            angles=(ManualAngle((1, 0, 2), 45.0, 1.57),),
        )
    assert exc.value.code == "incomplete_manual_bond_coverage"

    with pytest.raises(MetalAssignmentValidationError) as exc:
        _prepare(molecule, angles=())
    assert exc.value.code == "incomplete_manual_angle_coverage"


def test_manual_bonded_rejects_tampering_and_coordinate_drift():
    molecule = _load()
    plan = _prepare(molecule)
    tampered_reference = replace(
        plan.reference,
        bonds=(
            ManualBond((0, 1), 999.0, 2.1),
            ManualBond((0, 2), 140.0, 2.0),
        ),
    )
    tampered_plan = replace(plan, reference=tampered_reference)

    with pytest.raises(MetalAssignmentValidationError) as exc:
        apply_manual_bonded_assignment(molecule, tampered_plan)
    assert exc.value.code == "stale_manual_bonded_artifact_hash"

    molecule.atoms[1].x += 0.1
    with pytest.raises(MetalAssignmentValidationError) as exc:
        apply_manual_bonded_assignment(molecule, plan)
    assert exc.value.code == "manual_bonded_input_mismatch"


def test_manual_bonded_same_residue_coordination_uses_local_bonds():
    molecule = _load(SAME_RESIDUE_PDB)

    result = apply_manual_bonded_assignment(molecule, _prepare(molecule))

    assert result.molecule.coordination_bonds == [(0, 1), (0, 2)]
    assert result.molecule.explicit_bonds == []
    assert result.molecule.residue_links == []


def test_manual_bonded_parameters_reach_raw_and_bundle_savers(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _load()
    atom_parameters = (
        AtomParameterUpdate(0, "manual_Zn_export", 65.4, "fixture"),
        AtomParameterUpdate(1, "manual_N_export", 14.01, "fixture"),
        AtomParameterUpdate(2, "manual_O_export", 16.0, "fixture"),
    )
    lj_parameters = tuple(
        LJParameter(atom_type, atom_type, epsilon, rmin, "fixture")
        for atom_type, epsilon, rmin in (
            ("manual_Zn_export", 0.0125, 1.1),
            ("manual_N_export", 0.025, 1.5),
            ("manual_O_export", 0.03, 1.4),
        )
    )
    plan = prepare_manual_bonded_assignment(
        molecule,
        metal_sites=(MetalSite(0, "Zn", 2),),
        electronic_state=ElectronicState(0, 1),
        coordination_edges=((0, 1), (0, 2)),
        bonds=(
            ManualBond((0, 1), 120.0, 2.1),
            ManualBond((0, 2), 140.0, 2.0),
        ),
        angles=(ManualAngle((1, 0, 2), 45.0, 1.5707963267948966),),
        atom_parameters=atom_parameters,
        lj_parameters=lj_parameters,
    )
    applied = apply_manual_bonded_assignment(molecule, plan).molecule

    Xponge.Save_SPONGE_Input(applied, "manual", dirname=str(tmp_path / "raw"))
    Xponge.save_sponge_input_bundle(applied, "manual", tmp_path / "bundle")

    raw_bonds = (tmp_path / "raw" / "manual_bond.txt").read_text().splitlines()
    raw_angles = (tmp_path / "raw" / "manual_angle.txt").read_text().splitlines()
    assert raw_bonds[1:] == [
        "0 1 120.000000 2.100000",
        "0 2 140.000000 2.000000",
    ]
    assert raw_angles[1] == "1 0 2 45.000000 1.570796"
    with h5py.File(
        tmp_path / "bundle" / "manual_topology.spgt.h5", "r"
    ) as handle:
        assert handle["/forcefield/bond/k"][:].tolist() == pytest.approx(
            [120.0, 140.0]
        )
        assert handle["/forcefield/angle/k"][:].tolist() == pytest.approx([45.0])


def _read_mdout_first_frame(path: Path):
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise AssertionError(f"SPONGE did not write a frame to {path}")
    return {
        key: float(value)
        for key, value in zip(lines[0].split(), lines[1].split())
        if key != "step"
    }


@pytest.mark.skipif(
    not os.environ.get("SPONGE_EXECUTABLE"),
    reason="set SPONGE_EXECUTABLE to run the manual-bonded numerical gate",
)
def test_manual_bonded_raw_and_bundle_roundtrip_match_sponge(tmp_path):
    executable = Path(os.environ["SPONGE_EXECUTABLE"]).resolve()
    if not executable.is_file():
        pytest.fail(f"SPONGE_EXECUTABLE does not exist: {executable}")
    molecule = _load()
    molecule.set_box_padding(12.0)
    atom_parameters = (
        AtomParameterUpdate(0, "e2e_Zn", 65.4, "e2e"),
        AtomParameterUpdate(1, "e2e_N", 14.01, "e2e"),
        AtomParameterUpdate(2, "e2e_O", 16.0, "e2e"),
    )
    lj_parameters = (
        LJParameter("e2e_Zn", "e2e_Zn", 0.0125, 1.1, "e2e"),
        LJParameter("e2e_N", "e2e_N", 0.025, 1.5, "e2e"),
        LJParameter("e2e_O", "e2e_O", 0.03, 1.4, "e2e"),
    )
    plan = prepare_manual_bonded_assignment(
        molecule,
        metal_sites=(MetalSite(0, "Zn", 2),),
        electronic_state=ElectronicState(0, 1),
        coordination_edges=((0, 1), (0, 2)),
        bonds=(
            ManualBond((0, 1), 120.0, 2.1),
            ManualBond((0, 2), 140.0, 2.0),
        ),
        angles=(ManualAngle((1, 0, 2), 45.0, 1.5707963267948966),),
        atom_parameters=atom_parameters,
        lj_parameters=lj_parameters,
    )
    applied = apply_manual_bonded_assignment(molecule, plan).molecule
    raw_dir = tmp_path / "raw"
    bundle_dir = tmp_path / "bundle"
    roundtrip_dir = tmp_path / "roundtrip"
    Xponge.Save_SPONGE_Input(applied, "manual", dirname=str(raw_dir))
    Xponge.save_sponge_input_bundle(applied, "manual", bundle_dir)
    convert_bundle_to_legacy(bundle_dir, roundtrip_dir, prefix="manual")

    keys = (
        "residue",
        "resname",
        "atom_name",
        "atom_type_name",
        "mass",
        "charge",
        "coordinate",
        "LJ",
        "bond",
        "angle",
        "dihedral",
        "exclude",
        "nb14",
    )
    mdin = "\n".join(
        [
            'mode = "minimization"',
            "step_limit = 0",
            "cutoff = 8.0",
            "print_zeroth_frame = 1",
            "write_mdout_interval = 1",
            "write_information_interval = 1",
            *[f'{key}_in_file = "manual_{key}.txt"' for key in keys],
        ]
    )
    frames = []
    for run_dir in (raw_dir, roundtrip_dir):
        (run_dir / "mdin.spg.toml").write_text(mdin + "\n")
        completed = subprocess.run(
            [str(executable), "-mdin", "mdin.spg.toml"],
            cwd=run_dir,
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        assert completed.returncode == 0, (
            f"SPONGE failed in {run_dir}\nstdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
        frames.append(_read_mdout_first_frame(run_dir / "mdout.txt"))
    assert frames[0].keys() == frames[1].keys()
    for key in frames[0]:
        assert frames[1][key] == pytest.approx(
            frames[0][key], abs=1.0e-5, rel=1.0e-7
        ), key
