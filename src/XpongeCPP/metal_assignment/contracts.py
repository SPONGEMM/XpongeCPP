"""Immutable contracts for molecule-first metal assignment."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from typing import Any
from uuid import uuid4


SCHEMA_VERSION = 1


class MetalAssignmentValidationError(ValueError):
    """A stable validation error with a machine-readable code."""

    def __init__(self, code: str, message: str, *, path: str = ""):
        super().__init__(message)
        self.code = code
        self.path = path


def _canonical_hash(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class ElectronicState:
    """Explicit whole-model charge and spin multiplicity."""

    total_charge: int
    spin_multiplicity: int

    def validate(self) -> None:
        if self.spin_multiplicity < 1:
            raise MetalAssignmentValidationError(
                "invalid_spin_multiplicity",
                "spin_multiplicity must be positive",
                path="electronic_state.spin_multiplicity",
            )


@dataclass(frozen=True, slots=True)
class MetalSite:
    """One explicitly mapped metal atom."""

    atom_id: int
    element: str
    formal_charge: int


@dataclass(frozen=True, slots=True)
class ChargeUpdate:
    """One planned atom-charge replacement."""

    atom_id: int
    charge: float
    source: str = "constrained_resp"


@dataclass(frozen=True, slots=True)
class ChargeLedgerEntry:
    """Auditable charge delta for one parent atom."""

    atom_id: int
    old_charge: float
    new_charge: float
    delta: float
    source: str


@dataclass(frozen=True, slots=True)
class MetalAssignmentRequest:
    """Hash-closed, explicit molecule-first request."""

    metal_sites: tuple[MetalSite, ...]
    electronic_state: ElectronicState
    coordination_edges: tuple[tuple[int, int], ...]
    charge_updates: tuple[ChargeUpdate, ...] = ()
    interaction_model: str = "bonded"
    request_id: str = field(default_factory=lambda: str(uuid4()))
    schema_version: int = SCHEMA_VERSION

    def payload(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def request_hash(self) -> str:
        return _canonical_hash(self.payload())


@dataclass(frozen=True, slots=True)
class MetalAssignmentPlan:
    """Side-effect-free overlay plan bound to one molecule state."""

    request: MetalAssignmentRequest
    input_hash: str
    topology_hash: str
    charge_ledger: tuple[ChargeLedgerEntry, ...]
    link_overlay: tuple[tuple[int, int], ...]
    plan_hash: str = ""
    schema_version: int = SCHEMA_VERSION

    def payload(self) -> dict[str, Any]:
        data = asdict(self)
        data.pop("plan_hash", None)
        return data

    def computed_hash(self) -> str:
        return _canonical_hash(self.payload())

    def validate_hash(self) -> None:
        if not self.plan_hash or self.plan_hash != self.computed_hash():
            raise MetalAssignmentValidationError(
                "stale_plan_hash",
                "metal-assignment plan hash does not match its payload",
                path="plan_hash",
            )


@dataclass(frozen=True, slots=True)
class MetalAssignmentResult:
    """Applied molecule plus immutable provenance and audit facts."""

    molecule: Any
    plan: MetalAssignmentPlan
    result_input_hash: str
    result_topology_hash: str
    applied_charge_atom_ids: tuple[int, ...]
    applied_coordination_edges: tuple[tuple[int, int], ...]
    inplace: bool


def validate_finite_charge(update: ChargeUpdate, index: int) -> None:
    if not math.isfinite(float(update.charge)):
        raise MetalAssignmentValidationError(
            "nonfinite_charge",
            "metal-assignment charge updates must be finite",
            path=f"charge_updates[{index}].charge",
        )


__all__ = [
    "ChargeLedgerEntry",
    "ChargeUpdate",
    "ElectronicState",
    "MetalAssignmentPlan",
    "MetalAssignmentRequest",
    "MetalAssignmentResult",
    "MetalAssignmentValidationError",
    "MetalSite",
]
