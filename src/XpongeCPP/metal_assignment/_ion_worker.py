"""Isolated 12-6 ion-parameter worker."""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import importlib
import importlib.metadata
import json
import math
import re
import sys
from typing import Any, Mapping


WORKER_PROTOCOL_VERSION = 1
SUPPORTED_WATER_MODELS = frozenset({"tip3p", "tip4pew", "spce", "opc"})


class WorkerInputError(ValueError):
    """Invalid private ion worker request."""


def _canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    ).hexdigest()


def _strict_object(value: Any, required: set[str], path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise WorkerInputError(f"{path}: expected object")
    if set(value) != required:
        missing = required - set(value)
        unknown = set(value) - required
        detail = f"missing={','.join(sorted(missing))};unknown={','.join(sorted(unknown))}"
        raise WorkerInputError(f"{path}: {detail}")
    return value


def _provider_version() -> str:
    try:
        return importlib.metadata.version("Xponge")
    except importlib.metadata.PackageNotFoundError:
        return "source-tree"


def _normalize_element(value: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[A-Z][a-z]?", value):
        raise WorkerInputError(f"invalid element {value!r}")
    return value


def _atom_type_element(type_name: str) -> str | None:
    match = re.match(r"([A-Z][a-z]?)", type_name)
    return match.group(1) if match else None


def _execute(value: Any) -> dict[str, Any]:
    request = _strict_object(
        value,
        {"protocol_version", "topology_hash", "water_model", "metals"},
        "request",
    )
    if request["protocol_version"] != WORKER_PROTOCOL_VERSION:
        raise WorkerInputError("unsupported protocol_version")
    water_model = request["water_model"]
    if water_model not in SUPPORTED_WATER_MODELS:
        raise WorkerInputError(f"unsupported water_model {water_model!r}")
    if not isinstance(request["topology_hash"], str) or len(request["topology_hash"]) != 64:
        raise WorkerInputError("topology_hash: expected sha256")
    if not isinstance(request["metals"], list) or not request["metals"]:
        raise WorkerInputError("metals: expected non-empty array")
    metals: list[Mapping[str, Any]] = []
    metal_ids: set[str] = set()
    for index, raw_metal in enumerate(request["metals"]):
        metal = _strict_object(raw_metal, {"external_id", "element", "formal_charge"}, f"metals[{index}]")
        if not isinstance(metal["external_id"], str) or not metal["external_id"]:
            raise WorkerInputError(f"metals[{index}].external_id: expected non-empty string")
        if metal["external_id"] in metal_ids:
            raise WorkerInputError(f"metals[{index}].external_id: duplicate")
        metal_ids.add(metal["external_id"])
        _normalize_element(metal["element"])
        if (
            isinstance(metal["formal_charge"], bool)
            or not isinstance(metal["formal_charge"], int)
            or metal["formal_charge"] <= 0
        ):
            raise WorkerInputError(f"metals[{index}].formal_charge: expected positive integer")
        metals.append(metal)

    importlib.import_module(f"XpongeCPP.forcefield.amber.{water_model}")
    from XpongeCPP import has_template
    from XpongeCPP.helper import ResidueType
    from ._native_parameters import extract_native_parameters

    parameters: dict[str, Any] = {}
    source = f"xponge:ion-12-6:{water_model}:{_provider_version()}"
    for metal in metals:
        template_name = metal["element"].upper()
        if metal["formal_charge"] != 1:
            template_name += str(metal["formal_charge"])
        if not has_template(template_name):
            raise WorkerInputError(
                f"no unique ion template for {metal['element']}({metal['formal_charge']:+d}); "
                "candidate_count=0"
            )
        residue_type = ResidueType.get_type(template_name)
        if len(residue_type.atoms) != 1:
            raise WorkerInputError(
                f"ion template {template_name!r} is not monatomic"
            )
        atom = residue_type.atoms[0]
        element = _atom_type_element(str(atom.type))
        charge = float(atom.charge)
        if (
            element != metal["element"]
            or not math.isfinite(charge)
            or not math.isclose(
                charge, float(metal["formal_charge"]), abs_tol=1.0e-8
            )
        ):
            raise WorkerInputError(
                f"ion template {template_name!r} does not match requested state"
            )
        lj_parameters, bonded_records = extract_native_parameters(
            residue_type, [metal["external_id"]], source
        )
        if bonded_records:
            raise WorkerInputError(
                f"ion template {template_name!r} unexpectedly has bonded terms"
            )
        parameters[metal["external_id"]] = {
            "template_name": residue_type.name,
            "atom_type": str(atom.type),
            "charge": float(atom.charge),
            "mass": float(atom.mass),
            "lj": lj_parameters[metal["external_id"]],
        }
    response = {
        "ok": True,
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "topology_hash": request["topology_hash"],
        "water_model": water_model,
        "provider_version": _provider_version(),
        "parameters": parameters,
    }
    response["response_hash"] = _canonical_hash(response)
    return response


def main() -> int:
    try:
        value = json.load(sys.stdin)
        with redirect_stdout(sys.stderr):
            response = _execute(value)
        json.dump(response, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sys.stdout.write("\n")
        return 0
    except Exception as exc:
        json.dump(
            {
                "ok": False,
                "protocol_version": WORKER_PROTOCOL_VERSION,
                "error": {"type": type(exc).__name__, "message": str(exc)},
            },
            sys.stdout,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        sys.stdout.write("\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
