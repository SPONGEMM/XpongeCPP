from dataclasses import replace
from io import StringIO

import pytest
import XpongeCPP as Xponge
from XpongeCPP.metal_assignment import (
    AtomParameterUpdate,
    ElectronicState,
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


def test_manual_bonded_same_residue_coordination_uses_explicit_bonds():
    molecule = _load(SAME_RESIDUE_PDB)

    result = apply_manual_bonded_assignment(molecule, _prepare(molecule))

    assert result.molecule.explicit_bonds == [[0, 1], [0, 2]]
    assert result.molecule.residue_links == []


def test_manual_bonded_parameters_reach_raw_and_bundle_savers(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _load()
    atom_parameters = (
        AtomParameterUpdate(0, "manual_Zn_export", 65.4, "fixture"),
        AtomParameterUpdate(1, "manual_N_export", 14.01, "fixture"),
        AtomParameterUpdate(2, "manual_O_export", 16.0, "fixture"),
    )
    for atom_type, epsilon, rmin in (
        ("manual_Zn_export", 0.0125, 1.1),
        ("manual_N_export", 0.025, 1.5),
        ("manual_O_export", 0.03, 1.4),
    ):
        Xponge.register_amber_lj_parameter(
            atom_type, atom_type, epsilon, rmin
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
