"""Explicit-reference manual bonded metal-assignment mode."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math
from typing import Iterable

from .contracts import (
    AngleParameter,
    AtomParameterUpdate,
    BondParameter,
    ElectronicState,
    LJParameter,
    MetalAssignmentPlan,
    MetalAssignmentResult,
    MetalAssignmentValidationError,
    MetalSite,
    _canonical_hash,
)
from .molecule_api import (
    _normalize_edges,
    apply_metal_assignment,
    build_metal_parameter_overlay,
    molecule_input_hash,
    molecule_topology_hash,
    prepare_metal_assignment,
)


MANUAL_BONDED_SOURCE = "manual_bonded:explicit_reference_geometry"


@dataclass(frozen=True, slots=True)
class ManualBond:
    atom_ids: tuple[int, int]
    force_constant: float
    equilibrium_length: float


@dataclass(frozen=True, slots=True)
class ManualAngle:
    atom_ids: tuple[int, int, int]
    force_constant: float
    equilibrium_angle: float


@dataclass(frozen=True, slots=True)
class ManualBondedReference:
    """Hash-closed explicit geometry and force-constant artifact."""

    input_hash: str
    topology_hash: str
    bonds: tuple[ManualBond, ...]
    angles: tuple[ManualAngle, ...]
    equilibrium_geometry_source: str = "frozen_current_geometry"
    coordinate_unit: str = "angstrom"
    angle_unit: str = "radian"
    provider_revision: str = "explicit-reference-geometry-v1"
    artifact_hash: str = ""
    schema_version: int = 1

    def payload(self) -> dict:
        data = asdict(self)
        data.pop("artifact_hash", None)
        return data

    def computed_hash(self) -> str:
        return _canonical_hash(self.payload())

    def validate_hash(self) -> None:
        if not self.artifact_hash or self.artifact_hash != self.computed_hash():
            raise MetalAssignmentValidationError(
                "stale_manual_bonded_artifact_hash",
                "manual bonded artifact hash does not match its payload",
                path="reference.artifact_hash",
            )


@dataclass(frozen=True, slots=True)
class ManualBondedPlan:
    """Manual reference artifact bound to a generic assignment plan."""

    assignment_plan: MetalAssignmentPlan
    reference: ManualBondedReference
    plan_hash: str = ""
    schema_version: int = 1

    def payload(self) -> dict:
        return {
            "assignment_plan_hash": self.assignment_plan.plan_hash,
            "reference_artifact_hash": self.reference.artifact_hash,
            "schema_version": self.schema_version,
        }

    def computed_hash(self) -> str:
        return _canonical_hash(self.payload())

    def validate_hash(self) -> None:
        if not self.plan_hash or self.plan_hash != self.computed_hash():
            raise MetalAssignmentValidationError(
                "stale_manual_bonded_plan_hash",
                "manual bonded plan hash does not match its payload",
                path="plan_hash",
            )


def _canonical_angle(atom_ids: tuple[int, int, int]) -> tuple[int, int, int]:
    atom1, center, atom3 = map(int, atom_ids)
    return min(atom1, atom3), center, max(atom1, atom3)


def _expected_angles(
    metal_sites: tuple[MetalSite, ...], edges: tuple[tuple[int, int], ...]
) -> set[tuple[int, int, int]]:
    donors_by_metal = {int(site.atom_id): [] for site in metal_sites}
    for atom1, atom2 in edges:
        metal = atom1 if atom1 in donors_by_metal else atom2
        donor = atom2 if metal == atom1 else atom1
        donors_by_metal[metal].append(donor)
    expected = set()
    for metal, donors in donors_by_metal.items():
        donors = sorted(donors)
        for index, donor1 in enumerate(donors):
            for donor2 in donors[index + 1 :]:
                expected.add((donor1, metal, donor2))
    return expected


def _validate_manual_reference(
    molecule,
    reference: ManualBondedReference,
    *,
    metal_sites: tuple[MetalSite, ...],
    coordination_edges: tuple[tuple[int, int], ...],
) -> None:
    reference.validate_hash()
    if reference.input_hash != molecule_input_hash(molecule):
        raise MetalAssignmentValidationError(
            "manual_bonded_input_mismatch",
            "manual bonded artifact targets a different molecule input",
            path="reference.input_hash",
        )
    if reference.topology_hash != molecule_topology_hash(molecule):
        raise MetalAssignmentValidationError(
            "manual_bonded_topology_mismatch",
            "manual bonded artifact targets a different topology",
            path="reference.topology_hash",
        )
    if reference.equilibrium_geometry_source != "frozen_current_geometry":
        raise MetalAssignmentValidationError(
            "invalid_manual_geometry_source",
            "manual bonded mode requires frozen_current_geometry",
        )
    if reference.coordinate_unit != "angstrom" or reference.angle_unit != "radian":
        raise MetalAssignmentValidationError(
            "invalid_manual_units",
            "manual bonded mode requires angstrom lengths and radian angles",
        )
    bond_keys = []
    for index, term in enumerate(reference.bonds):
        if len(term.atom_ids) != 2:
            raise MetalAssignmentValidationError(
                "invalid_manual_bond_atoms",
                "manual bond requires two atom ids",
                path=f"reference.bonds[{index}].atom_ids",
            )
        key = tuple(sorted(map(int, term.atom_ids)))
        if (
            key[0] == key[1]
            or not math.isfinite(float(term.force_constant))
            or not 0.0 < float(term.force_constant) <= 10000.0
            or not math.isfinite(float(term.equilibrium_length))
            or float(term.equilibrium_length) <= 0.0
        ):
            raise MetalAssignmentValidationError(
                "invalid_manual_bond",
                "manual bond requires positive finite constants and length",
                path=f"reference.bonds[{index}]",
            )
        bond_keys.append(key)
    expected_bonds = set(coordination_edges)
    if set(bond_keys) != expected_bonds or len(bond_keys) != len(set(bond_keys)):
        raise MetalAssignmentValidationError(
            "incomplete_manual_bond_coverage",
            "manual bonds must cover every coordination edge exactly once",
            path="reference.bonds",
        )
    angle_keys = []
    for index, term in enumerate(reference.angles):
        if len(term.atom_ids) != 3:
            raise MetalAssignmentValidationError(
                "invalid_manual_angle_atoms",
                "manual angle requires three atom ids",
                path=f"reference.angles[{index}].atom_ids",
            )
        key = _canonical_angle(term.atom_ids)
        if (
            len(set(key)) != 3
            or not math.isfinite(float(term.force_constant))
            or not 0.0 < float(term.force_constant) <= 10000.0
            or not math.isfinite(float(term.equilibrium_angle))
            or not 0.0 < float(term.equilibrium_angle) < math.pi
        ):
            raise MetalAssignmentValidationError(
                "invalid_manual_angle",
                "manual angle requires positive finite constants and angle in (0, pi)",
                path=f"reference.angles[{index}]",
            )
        angle_keys.append(key)
    expected_angles = _expected_angles(metal_sites, coordination_edges)
    if set(angle_keys) != expected_angles or len(angle_keys) != len(set(angle_keys)):
        raise MetalAssignmentValidationError(
            "incomplete_manual_angle_coverage",
            "manual angles must cover every donor pair around each metal exactly once",
            path="reference.angles",
        )


def prepare_manual_bonded_assignment(
    molecule,
    *,
    metal_sites: Iterable[MetalSite],
    electronic_state: ElectronicState,
    coordination_edges: Iterable[tuple[int, int]],
    bonds: Iterable[ManualBond],
    angles: Iterable[ManualAngle] = (),
    atom_parameters: Iterable[AtomParameterUpdate] = (),
    lj_parameters: Iterable[LJParameter] = (),
) -> ManualBondedPlan:
    """Prepare manual bonded terms without QM, RESP, or parent mutation."""

    sites = tuple(metal_sites)
    edges = _normalize_edges(coordination_edges)
    reference = ManualBondedReference(
        input_hash=molecule_input_hash(molecule),
        topology_hash=molecule_topology_hash(molecule),
        bonds=tuple(bonds),
        angles=tuple(angles),
    )
    reference = replace(reference, artifact_hash=reference.computed_hash())
    _validate_manual_reference(
        molecule,
        reference,
        metal_sites=sites,
        coordination_edges=edges,
    )
    overlay = build_metal_parameter_overlay(
        molecule,
        atom_parameters=tuple(atom_parameters),
        bond_parameters=tuple(
            BondParameter(
                atom_ids=tuple(map(int, term.atom_ids)),
                force_constant=float(term.force_constant),
                equilibrium_length=float(term.equilibrium_length),
                source=MANUAL_BONDED_SOURCE,
            )
            for term in reference.bonds
        ),
        angle_parameters=tuple(
            AngleParameter(
                atom_ids=tuple(map(int, term.atom_ids)),
                force_constant=float(term.force_constant),
                equilibrium_angle=float(term.equilibrium_angle),
                source=MANUAL_BONDED_SOURCE,
            )
            for term in reference.angles
        ),
        lj_parameters=tuple(lj_parameters),
        parameter_source=MANUAL_BONDED_SOURCE,
    )
    assignment_plan = prepare_metal_assignment(
        molecule,
        metal_sites=sites,
        electronic_state=electronic_state,
        coordination_edges=edges,
        interaction_model="bonded",
        parameter_overlay=overlay,
    )
    plan = ManualBondedPlan(
        assignment_plan=assignment_plan,
        reference=reference,
    )
    return replace(plan, plan_hash=plan.computed_hash())


def apply_manual_bonded_assignment(
    molecule, plan: ManualBondedPlan, *, inplace: bool = False
) -> MetalAssignmentResult:
    """Validate the explicit reference and transactionally apply its overlay."""

    plan.validate_hash()
    plan.assignment_plan.validate_hash()
    request = plan.assignment_plan.request
    _validate_manual_reference(
        molecule,
        plan.reference,
        metal_sites=request.metal_sites,
        coordination_edges=request.coordination_edges,
    )
    return apply_metal_assignment(
        molecule,
        plan.assignment_plan,
        inplace=inplace,
    )


__all__ = [
    "MANUAL_BONDED_SOURCE",
    "ManualAngle",
    "ManualBond",
    "ManualBondedPlan",
    "ManualBondedReference",
    "apply_manual_bonded_assignment",
    "prepare_manual_bonded_assignment",
]
