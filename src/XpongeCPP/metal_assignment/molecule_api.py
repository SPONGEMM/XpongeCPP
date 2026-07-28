"""Molecule-first planning and transactional application."""

from __future__ import annotations

from dataclasses import replace
import math
from typing import Iterable

from .contracts import (
    AngleParameter,
    AtomParameterUpdate,
    BondParameter,
    ChargeLedgerEntry,
    ChargeUpdate,
    ElectronicState,
    MetalAssignmentPlan,
    MetalAssignmentRequest,
    MetalAssignmentResult,
    MetalAssignmentValidationError,
    MetalParameterOverlay,
    MetalSite,
    _canonical_hash,
    validate_finite_charge,
)


def _atom_residue_ids(molecule) -> tuple[int, ...]:
    residue_ids = [-1] * int(molecule.atom_count)
    for residue in molecule.residues:
        for atom in residue.atoms:
            residue_ids[int(atom.index)] = int(residue.index)
    if any(residue_id < 0 for residue_id in residue_ids):
        raise MetalAssignmentValidationError(
            "invalid_parent_molecule",
            "parent molecule contains an atom without a residue",
        )
    return tuple(residue_ids)


def molecule_topology_hash(molecule) -> str:
    """Hash stable topology identity without coordinates or charges."""

    payload = {
        "atom_count": int(molecule.atom_count),
        "residue_count": int(molecule.residue_count),
        "atom_residue_ids": _atom_residue_ids(molecule),
        "residues": tuple(
            (
                int(residue.index),
                str(residue.name),
                int(residue.atom_count),
            )
            for residue in molecule.residues
        ),
        "explicit_bonds": tuple(
            sorted(tuple(sorted(map(int, edge))) for edge in molecule.explicit_bonds)
        ),
        "residue_links": tuple(
            sorted(tuple(sorted(map(int, edge))) for edge in molecule.residue_links)
        ),
    }
    return _canonical_hash(payload)


def molecule_input_hash(molecule) -> str:
    """Hash topology plus mutable atom state consumed by a plan."""

    payload = {
        "topology_hash": molecule_topology_hash(molecule),
        "atoms": tuple(
            (
                int(atom.index),
                str(atom.element),
                str(atom.name),
                str(atom.type),
                float(atom.mass),
                float(atom.charge),
                float(atom.x),
                float(atom.y),
                float(atom.z),
            )
            for atom in molecule.atoms
        ),
        "box": {
            "length": tuple(map(float, molecule.box_length)),
            "origin": tuple(map(float, molecule.box_origin)),
            "angle": tuple(map(float, molecule.box_angle)),
        },
    }
    return _canonical_hash(payload)


def _normalize_edges(edges: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    normalized = []
    for index, edge in enumerate(edges):
        if len(edge) != 2:
            raise MetalAssignmentValidationError(
                "invalid_coordination_edge",
                "coordination edges must contain exactly two atom ids",
                path=f"coordination_edges[{index}]",
            )
        atom1, atom2 = map(int, edge)
        if atom1 == atom2:
            raise MetalAssignmentValidationError(
                "self_coordination_edge",
                "coordination edges cannot be self edges",
                path=f"coordination_edges[{index}]",
            )
        normalized.append(tuple(sorted((atom1, atom2))))
    if len(normalized) != len(set(normalized)):
        raise MetalAssignmentValidationError(
            "duplicate_coordination_edge",
            "coordination edges must be unique",
            path="coordination_edges",
        )
    return tuple(sorted(normalized))


def _validate_parameter_overlay(
    molecule, overlay: MetalParameterOverlay, topology_hash: str
) -> None:
    overlay.validate_hash()
    if overlay.topology_hash != topology_hash:
        raise MetalAssignmentValidationError(
            "parameter_overlay_topology_mismatch",
            "metal parameter overlay targets a different topology",
            path="parameter_overlay.topology_hash",
        )
    if not overlay.parameter_source:
        raise MetalAssignmentValidationError(
            "missing_parameter_source",
            "metal parameter overlay requires a parameter source",
            path="parameter_overlay.parameter_source",
        )
    if overlay.precedence <= 0:
        raise MetalAssignmentValidationError(
            "invalid_parameter_precedence",
            "metal parameter overlay precedence must be positive",
            path="parameter_overlay.precedence",
        )
    atom_count = int(molecule.atom_count)
    atom_parameter_ids = []
    for index, update in enumerate(overlay.atom_parameters):
        atom_id = int(update.atom_id)
        if atom_id < 0 or atom_id >= atom_count:
            raise MetalAssignmentValidationError(
                "parameter_atom_out_of_range",
                "atom parameter update references an atom outside the molecule",
                path=f"parameter_overlay.atom_parameters[{index}].atom_id",
            )
        atom_parameter_ids.append(atom_id)
        if update.atom_type is None and update.mass is None:
            raise MetalAssignmentValidationError(
                "empty_atom_parameter_update",
                "atom parameter update must set atom_type or mass",
                path=f"parameter_overlay.atom_parameters[{index}]",
            )
        if update.atom_type is not None and not str(update.atom_type).strip():
            raise MetalAssignmentValidationError(
                "invalid_atom_type",
                "atom parameter type must not be empty",
                path=f"parameter_overlay.atom_parameters[{index}].atom_type",
            )
        if update.mass is not None and (
            not math.isfinite(float(update.mass)) or float(update.mass) <= 0.0
        ):
            raise MetalAssignmentValidationError(
                "invalid_atom_mass",
                "atom parameter mass must be finite and positive",
                path=f"parameter_overlay.atom_parameters[{index}].mass",
            )
        if not str(update.source):
            raise MetalAssignmentValidationError(
                "missing_atom_parameter_source",
                "atom parameter update requires a source",
                path=f"parameter_overlay.atom_parameters[{index}].source",
            )
    if len(atom_parameter_ids) != len(set(atom_parameter_ids)):
        raise MetalAssignmentValidationError(
            "duplicate_atom_parameter_update",
            "atom parameter update ids must be unique",
            path="parameter_overlay.atom_parameters",
        )
    bond_keys = []
    for index, term in enumerate(overlay.bond_parameters):
        if len(term.atom_ids) != 2:
            raise MetalAssignmentValidationError(
                "invalid_bond_parameter_atoms",
                "bond parameter requires exactly two atom ids",
                path=f"parameter_overlay.bond_parameters[{index}].atom_ids",
            )
        atom1, atom2 = map(int, term.atom_ids)
        key = tuple(sorted((atom1, atom2)))
        if atom1 == atom2 or key[0] < 0 or key[1] >= atom_count:
            raise MetalAssignmentValidationError(
                "invalid_bond_parameter_atoms",
                "bond parameter atom ids must be distinct and in range",
                path=f"parameter_overlay.bond_parameters[{index}].atom_ids",
            )
        if (
            not math.isfinite(float(term.force_constant))
            or float(term.force_constant) < 0.0
            or not math.isfinite(float(term.equilibrium_length))
            or float(term.equilibrium_length) <= 0.0
        ):
            raise MetalAssignmentValidationError(
                "invalid_bond_parameter",
                "bond force constant must be non-negative and length positive",
                path=f"parameter_overlay.bond_parameters[{index}]",
            )
        if not str(term.source):
            raise MetalAssignmentValidationError(
                "missing_bond_parameter_source",
                "bond parameter requires a source",
                path=f"parameter_overlay.bond_parameters[{index}].source",
            )
        bond_keys.append(key)
    if len(bond_keys) != len(set(bond_keys)):
        raise MetalAssignmentValidationError(
            "duplicate_bond_parameter",
            "bond parameter atom pairs must be unique",
            path="parameter_overlay.bond_parameters",
        )
    angle_keys = []
    for index, term in enumerate(overlay.angle_parameters):
        if len(term.atom_ids) != 3:
            raise MetalAssignmentValidationError(
                "invalid_angle_parameter_atoms",
                "angle parameter requires exactly three atom ids",
                path=f"parameter_overlay.angle_parameters[{index}].atom_ids",
            )
        atom1, center, atom3 = map(int, term.atom_ids)
        key = (min(atom1, atom3), center, max(atom1, atom3))
        if (
            len({atom1, center, atom3}) != 3
            or min(key) < 0
            or max(key) >= atom_count
        ):
            raise MetalAssignmentValidationError(
                "invalid_angle_parameter_atoms",
                "angle parameter atom ids must be distinct and in range",
                path=f"parameter_overlay.angle_parameters[{index}].atom_ids",
            )
        if (
            not math.isfinite(float(term.force_constant))
            or float(term.force_constant) < 0.0
            or not math.isfinite(float(term.equilibrium_angle))
            or not 0.0 < float(term.equilibrium_angle) < math.pi
        ):
            raise MetalAssignmentValidationError(
                "invalid_angle_parameter",
                "angle force constant must be non-negative and angle in (0, pi)",
                path=f"parameter_overlay.angle_parameters[{index}]",
            )
        if not str(term.source):
            raise MetalAssignmentValidationError(
                "missing_angle_parameter_source",
                "angle parameter requires a source",
                path=f"parameter_overlay.angle_parameters[{index}].source",
            )
        angle_keys.append(key)
    if len(angle_keys) != len(set(angle_keys)):
        raise MetalAssignmentValidationError(
            "duplicate_angle_parameter",
            "angle parameter atom triples must be unique",
            path="parameter_overlay.angle_parameters",
        )


def build_metal_parameter_overlay(
    molecule,
    *,
    atom_parameters: Iterable[AtomParameterUpdate] = (),
    bond_parameters: Iterable[BondParameter] = (),
    angle_parameters: Iterable[AngleParameter] = (),
    parameter_source: str = "metal_assignment",
    precedence: int = 100,
) -> MetalParameterOverlay:
    """Build and validate a hash-closed, molecule-local parameter overlay."""

    overlay = MetalParameterOverlay(
        topology_hash=molecule_topology_hash(molecule),
        atom_parameters=tuple(atom_parameters),
        bond_parameters=tuple(bond_parameters),
        angle_parameters=tuple(angle_parameters),
        parameter_source=str(parameter_source),
        precedence=int(precedence),
    )
    overlay = replace(overlay, overlay_hash=overlay.computed_hash())
    _validate_parameter_overlay(molecule, overlay, overlay.topology_hash)
    return overlay


def prepare_metal_assignment(
    molecule,
    *,
    metal_sites: Iterable[MetalSite],
    electronic_state: ElectronicState,
    coordination_edges: Iterable[tuple[int, int]],
    charge_updates: Iterable[ChargeUpdate] = (),
    interaction_model: str = "bonded",
    parameter_overlay: MetalParameterOverlay | None = None,
) -> MetalAssignmentPlan:
    """Create a validated plan without mutating the parent molecule."""

    if not molecule.validate():
        raise MetalAssignmentValidationError(
            "invalid_parent_molecule", "parent molecule is invalid"
        )
    sites = tuple(metal_sites)
    if not sites:
        raise MetalAssignmentValidationError(
            "missing_metal_sites", "at least one explicit metal site is required"
        )
    electronic_state.validate()
    if interaction_model not in {"bonded", "nonbonded_12_6"}:
        raise MetalAssignmentValidationError(
            "invalid_interaction_model",
            "interaction_model must be 'bonded' or 'nonbonded_12_6'",
            path="interaction_model",
        )
    atom_count = int(molecule.atom_count)
    site_ids = [int(site.atom_id) for site in sites]
    if len(site_ids) != len(set(site_ids)):
        raise MetalAssignmentValidationError(
            "duplicate_metal_site", "metal atom ids must be unique"
        )
    for index, site in enumerate(sites):
        if site.atom_id < 0 or site.atom_id >= atom_count:
            raise MetalAssignmentValidationError(
                "metal_atom_out_of_range",
                "metal atom id is outside the parent molecule",
                path=f"metal_sites[{index}].atom_id",
            )
        actual = str(molecule.atoms[site.atom_id].element).strip().lower()
        if actual != str(site.element).strip().lower():
            raise MetalAssignmentValidationError(
                "metal_element_mismatch",
                f"metal site element {site.element!r} does not match {actual!r}",
                path=f"metal_sites[{index}].element",
            )
    edges = _normalize_edges(coordination_edges)
    if interaction_model == "bonded" and not edges:
        raise MetalAssignmentValidationError(
            "missing_coordination_edges",
            "bonded metal assignment requires explicit coordination edges",
        )
    site_set = set(site_ids)
    residue_ids = _atom_residue_ids(molecule)
    for index, (atom1, atom2) in enumerate(edges):
        if atom1 < 0 or atom2 >= atom_count:
            raise MetalAssignmentValidationError(
                "coordination_atom_out_of_range",
                "coordination edge references an atom outside the molecule",
                path=f"coordination_edges[{index}]",
            )
        if (atom1 in site_set) == (atom2 in site_set):
            raise MetalAssignmentValidationError(
                "invalid_metal_coordination_edge",
                "each coordination edge must contain exactly one metal site",
                path=f"coordination_edges[{index}]",
            )
    updates = tuple(charge_updates)
    update_ids = [int(update.atom_id) for update in updates]
    if len(update_ids) != len(set(update_ids)):
        raise MetalAssignmentValidationError(
            "duplicate_charge_update", "charge update atom ids must be unique"
        )
    ledger = []
    for index, update in enumerate(updates):
        validate_finite_charge(update, index)
        if update.atom_id < 0 or update.atom_id >= atom_count:
            raise MetalAssignmentValidationError(
                "charge_atom_out_of_range",
                "charge update references an atom outside the molecule",
                path=f"charge_updates[{index}].atom_id",
            )
        old_charge = float(molecule.atoms[update.atom_id].charge)
        ledger.append(
            ChargeLedgerEntry(
                atom_id=int(update.atom_id),
                old_charge=old_charge,
                new_charge=float(update.charge),
                delta=float(update.charge) - old_charge,
                source=str(update.source),
            )
        )
    topology_hash = molecule_topology_hash(molecule)
    if parameter_overlay is not None:
        _validate_parameter_overlay(molecule, parameter_overlay, topology_hash)
    request = MetalAssignmentRequest(
        metal_sites=sites,
        electronic_state=electronic_state,
        coordination_edges=edges,
        charge_updates=updates,
        interaction_model=interaction_model,
    )
    plan = MetalAssignmentPlan(
        request=request,
        input_hash=molecule_input_hash(molecule),
        topology_hash=topology_hash,
        charge_ledger=tuple(ledger),
        link_overlay=edges,
        parameter_overlay=parameter_overlay,
    )
    return replace(plan, plan_hash=plan.computed_hash())


def apply_metal_assignment(
    molecule, plan: MetalAssignmentPlan, *, inplace: bool = False
) -> MetalAssignmentResult:
    """Apply a validated plan on a copy, committing only after full success."""

    plan.validate_hash()
    current_hash = molecule_input_hash(molecule)
    if current_hash != plan.input_hash:
        raise MetalAssignmentValidationError(
            "stale_parent_molecule",
            "parent molecule changed after the metal-assignment plan was made",
            path="input_hash",
        )
    if plan.link_overlay and molecule.has_topology_override:
        raise MetalAssignmentValidationError(
            "topology_override_conflict",
            "cannot apply coordination links while topology_override is active",
        )
    working = molecule.deepcopy()
    for entry in plan.charge_ledger:
        if not math.isclose(
            float(working.atoms[entry.atom_id].charge),
            entry.old_charge,
            abs_tol=1e-12,
        ):
            raise MetalAssignmentValidationError(
                "stale_charge_ledger",
                f"atom {entry.atom_id} charge no longer matches the plan",
            )
        working.atoms[entry.atom_id].charge = entry.new_charge
    overlay = plan.parameter_overlay
    if overlay is not None:
        _validate_parameter_overlay(working, overlay, plan.topology_hash)
        for update in overlay.atom_parameters:
            atom = working.atoms[update.atom_id]
            if update.atom_type is not None:
                atom.type = str(update.atom_type)
            if update.mass is not None:
                atom.mass = float(update.mass)
        for term in overlay.bond_parameters:
            working._set_bond_parameter_override(
                *term.atom_ids,
                float(term.force_constant),
                float(term.equilibrium_length),
                str(term.source),
            )
        for term in overlay.angle_parameters:
            working._set_angle_parameter_override(
                *term.atom_ids,
                float(term.force_constant),
                float(term.equilibrium_angle),
                str(term.source),
            )
    residue_ids = _atom_residue_ids(working)
    for atom1, atom2 in plan.link_overlay:
        if residue_ids[atom1] == residue_ids[atom2]:
            working.add_explicit_bond(atom1, atom2)
        else:
            working.add_residue_link(atom1, atom2)
    if not working.validate():
        raise MetalAssignmentValidationError(
            "invalid_applied_molecule",
            "metal-assignment overlay produced an invalid molecule",
        )
    result_input_hash = molecule_input_hash(working)
    result_topology_hash = molecule_topology_hash(working)
    if inplace:
        molecule._replace_from(working)
        published = molecule
    else:
        published = working
    return MetalAssignmentResult(
        molecule=published,
        plan=plan,
        result_input_hash=result_input_hash,
        result_topology_hash=result_topology_hash,
        applied_charge_atom_ids=tuple(
            entry.atom_id for entry in plan.charge_ledger
        ),
        applied_coordination_edges=plan.link_overlay,
        inplace=bool(inplace),
    )


__all__ = [
    "apply_metal_assignment",
    "build_metal_parameter_overlay",
    "molecule_input_hash",
    "molecule_topology_hash",
    "prepare_metal_assignment",
]
