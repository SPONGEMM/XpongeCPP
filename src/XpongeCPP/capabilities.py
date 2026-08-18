"""Machine-readable Xponge compatibility capabilities.

Only capabilities declared here are part of the stable compatibility contract.
An unlisted capability is unsupported until it is implemented and covered by
parity tests.  ``parity`` distinguishes an available implementation from one
whose numerical equivalence has completed the release gate.
"""

from __future__ import annotations

from copy import deepcopy


_CAPABILITIES = {
    "compat.import_xponge": {"status": "supported", "parity": "verified"},
    "forcefield.amber.protein": {"status": "supported", "parity": "verified"},
    "forcefield.amber.water_ion": {"status": "supported", "parity": "verified"},
    "forcefield.amber.gaff": {"status": "supported", "parity": "verified"},
    "forcefield.amber.gaff2": {"status": "supported", "parity": "verified"},
    "forcefield.amber.glycam": {"status": "supported", "parity": "verified"},
    "forcefield.amber.lipid17": {"status": "supported", "parity": "verified"},
    "forcefield.amber.lipid21": {"status": "supported", "parity": "verified"},
    "io.pdb": {"status": "supported", "parity": "verified"},
    "io.mmcif": {"status": "supported", "parity": "verified"},
    "io.mol2": {"status": "supported", "parity": "verified"},
    "io.sponge.raw": {"status": "supported", "parity": "verified"},
    "io.sponge.bundle": {"status": "supported", "parity": "verified"},
    "assignment.resp": {"status": "supported", "parity": "verified"},
    "metal_assignment.local_patch": {"status": "supported", "parity": "verified"},
    "fep.dual_topology": {"status": "supported", "parity": "verified"},
}


def capability_manifest():
    """Return an isolated JSON-serializable compatibility manifest."""

    return {
        "schema_version": 1,
        "implementation": "xpongecpp",
        "unlisted_status": "unsupported",
        "capabilities": deepcopy(_CAPABILITIES),
    }


def capability_status(capability_id):
    """Return the declared status for *capability_id*."""

    entry = _CAPABILITIES.get(str(capability_id))
    return "unsupported" if entry is None else entry["status"]


def require_capability(capability_id):
    """Fail explicitly when a requested compatibility capability is absent."""

    capability_id = str(capability_id)
    if capability_status(capability_id) != "supported":
        raise NotImplementedError(
            f"XpongeCPP compatibility capability is not supported: {capability_id}"
        )


__all__ = [
    "capability_manifest",
    "capability_status",
    "require_capability",
]
