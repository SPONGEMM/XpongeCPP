"""Molecule-first metal-assignment facade."""

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
)
from .molecule_api import (
    apply_metal_assignment,
    build_metal_parameter_overlay,
    molecule_input_hash,
    molecule_topology_hash,
    prepare_metal_assignment,
)
from .manual_bonded import (
    MANUAL_BONDED_SOURCE,
    ManualAngle,
    ManualBond,
    ManualBondedPlan,
    ManualBondedReference,
    apply_manual_bonded_assignment,
    prepare_manual_bonded_assignment,
)

__all__ = [
    "AngleParameter",
    "AtomParameterUpdate",
    "BondParameter",
    "ChargeLedgerEntry",
    "ChargeUpdate",
    "ElectronicState",
    "MetalAssignmentPlan",
    "MetalAssignmentRequest",
    "MetalAssignmentResult",
    "MetalAssignmentValidationError",
    "MetalParameterOverlay",
    "MetalSite",
    "MANUAL_BONDED_SOURCE",
    "ManualAngle",
    "ManualBond",
    "ManualBondedPlan",
    "ManualBondedReference",
    "apply_metal_assignment",
    "apply_manual_bonded_assignment",
    "build_metal_parameter_overlay",
    "molecule_input_hash",
    "molecule_topology_hash",
    "prepare_metal_assignment",
    "prepare_manual_bonded_assignment",
]
