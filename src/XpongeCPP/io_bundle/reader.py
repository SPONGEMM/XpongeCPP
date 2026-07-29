"""Validated HDF5 reader for bundled SPONGE input artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .case import BundleCase
from .errors import BundlePathError, BundleSchemaError, BundleValidationError


_INPUT_SCHEMA_VERSION = "sponge.input.v2"
_OUTPUT_SCHEMA_VERSION = "sponge.output.v2"
_SCHEMA_BINDINGS = {
    "topology.spgt.h5": ("/schema/name", "/schema/version", "sponge.topology.h5"),
    "protocol.spgp.h5": ("/schema/name", "/schema/version", "sponge.protocol.h5"),
    "restart.spgr.h5": ("/schema/name", "/schema/version", "sponge.restart.h5"),
    "trajectory.spg.h5md": (
        "/parameters/sponge/schema/name",
        "/parameters/sponge/schema/version",
        "sponge.output.h5md",
    ),
}


class BundleReader:
    """Keep bundle handles open while validating and reading typed datasets."""

    def __init__(self, case: BundleCase, *, strict: bool = True):
        self.case = case
        self.strict = strict
        self._handles: dict[str, Any] = {}
        self.warnings: list[str] = []

    def __enter__(self) -> "BundleReader":
        try:
            import h5py
        except ImportError as exc:  # pragma: no cover - required dependency
            raise BundleValidationError(
                "h5py is required to read bundled inputs"
            ) from exc

        for bundle_file in _SCHEMA_BINDINGS:
            path = self.case.path_for_bundle_file(bundle_file)
            if path is not None and path.is_file():
                self._handles[bundle_file] = h5py.File(path, "r")
        try:
            self.validate()
        except Exception:
            self.close()
            raise
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def close(self) -> None:
        for handle in self._handles.values():
            handle.close()
        self._handles.clear()

    def has_bundle_file(self, bundle_file: str) -> bool:
        return bundle_file in self._handles

    def contains(self, bundle_file: str, dataset_path: str) -> bool:
        handle = self._handles.get(bundle_file)
        return handle is not None and dataset_path in handle

    def read(self, bundle_file: str, dataset_path: str):
        handle = self._require_handle(bundle_file)
        if dataset_path not in handle:
            raise BundleValidationError(f"{bundle_file} is missing {dataset_path}")
        return np.asarray(handle[dataset_path][...])

    def read_scalar(self, bundle_file: str, dataset_path: str):
        value = self.read(bundle_file, dataset_path)
        return value.item() if value.shape == () else value.reshape(-1)[0].item()

    def read_text(self, bundle_file: str, dataset_path: str) -> str:
        handle = self._require_handle(bundle_file)
        if dataset_path not in handle:
            raise BundleValidationError(f"{bundle_file} is missing {dataset_path}")
        return _decode_h5_text(handle[dataset_path][()])

    def read_attribute(
        self, bundle_file: str, dataset_path: str, attribute: str
    ):
        """Read one HDF5 attribute from a validated artifact."""

        handle = self._require_handle(bundle_file)
        if dataset_path not in handle:
            raise BundleValidationError(f"{bundle_file} is missing {dataset_path}")
        if attribute not in handle[dataset_path].attrs:
            raise BundleValidationError(
                f"{bundle_file}:{dataset_path} is missing attribute {attribute!r}"
            )
        return _decode_h5_text(handle[dataset_path].attrs[attribute])

    def read_legacy_sidecars(self, bundle_file: str) -> dict[str, Path]:
        key_path = "/parameters/sponge/files/legacy_sidecars/key"
        value_path = "/parameters/sponge/files/legacy_sidecars/path"
        if not self.contains(bundle_file, key_path) and not self.contains(
            bundle_file, value_path
        ):
            return {}
        if not self.contains(bundle_file, key_path) or not self.contains(
            bundle_file, value_path
        ):
            raise BundleValidationError(
                f"{bundle_file} has an incomplete legacy sidecar table"
            )

        keys = [
            _decode_h5_text(value)
            for value in self.read(bundle_file, key_path).reshape(-1)
        ]
        paths = [
            _decode_h5_text(value)
            for value in self.read(bundle_file, value_path).reshape(-1)
        ]
        if len(keys) != len(paths):
            raise BundleValidationError(
                f"{bundle_file} legacy sidecar key/path lengths differ"
            )

        sidecars: dict[str, Path] = {}
        root = self.case.root.resolve()
        for key, raw_path in zip(keys, paths):
            path = (root / raw_path).resolve()
            try:
                path.relative_to(root)
            except ValueError as exc:
                raise BundlePathError(
                    f"legacy sidecar for {key} escapes bundle root: {raw_path}"
                ) from exc
            if not path.is_file():
                self._validation_problem(
                    f"legacy sidecar for {key} does not exist: {path}"
                )
                continue
            sidecars[key] = path
        return sidecars

    def validate(self) -> list[str]:
        for bundle_file, (
            name_path,
            version_path,
            expected_name,
        ) in _SCHEMA_BINDINGS.items():
            if bundle_file not in self._handles:
                continue
            self._validate_schema(
                bundle_file, name_path, version_path, expected_name
            )
            self.read_legacy_sidecars(bundle_file)
        self._validate_identity_uuid()
        self._validate_lineage()
        self._validate_atom_dimensions()
        self._validate_finalized_restart()
        self._validate_restart_load_policy()
        return list(self.warnings)

    def _validate_schema(
        self,
        bundle_file: str,
        name_path: str,
        version_path: str,
        expected_name: str,
    ) -> None:
        if not self.contains(bundle_file, name_path) or not self.contains(
            bundle_file, version_path
        ):
            self._schema_problem(f"{bundle_file} is missing schema metadata")
            return
        name = self.read_text(bundle_file, name_path)
        version = self.read_text(bundle_file, version_path)
        if name != expected_name:
            self._schema_problem(
                f"{bundle_file} schema name is {name!r}, "
                f"expected {expected_name!r}"
            )
        expected_version = (
            _OUTPUT_SCHEMA_VERSION
            if bundle_file == "trajectory.spg.h5md"
            else _INPUT_SCHEMA_VERSION
        )
        if version != expected_version:
            self._schema_problem(
                f"{bundle_file} schema version is {version!r}, "
                f"expected {expected_version!r}"
            )

    def _validate_identity_uuid(self) -> None:
        identities = {
            bundle_file: self.read_text(bundle_file, "/identity/uuid")
            for bundle_file in (
                "topology.spgt.h5",
                "protocol.spgp.h5",
                "restart.spgr.h5",
            )
            if self.contains(bundle_file, "/identity/uuid")
        }
        if len(set(identities.values())) > 1:
            raise BundleValidationError(
                f"bundle identity UUIDs do not match: {identities}"
            )

    def _validate_lineage(self) -> None:
        comparisons = (
            (
                "topology.spgt.h5",
                "/topology/topology_hash",
                "protocol.spgp.h5",
                "/protocol/topology_compatibility/topology_hash",
                "protocol topology compatibility hash",
            ),
            (
                "topology.spgt.h5",
                "/topology/topology_hash",
                "restart.spgr.h5",
                "/run/topology_hash",
                "restart topology hash",
            ),
            (
                "topology.spgt.h5",
                "/topology/atom_order_hash",
                "restart.spgr.h5",
                "/run/atom_order_hash",
                "restart atom-order hash",
            ),
            (
                "protocol.spgp.h5",
                "/identity/content_hash",
                "restart.spgr.h5",
                "/run/producer_protocol_hash",
                "restart producer protocol hash",
            ),
        )
        for left_file, left_path, right_file, right_path, label in comparisons:
            if not self.contains(left_file, left_path) or not self.contains(
                right_file, right_path
            ):
                continue
            if self.read_text(left_file, left_path) != self.read_text(
                right_file, right_path
            ):
                raise BundleValidationError(f"{label} does not match")

    def _validate_atom_dimensions(self) -> None:
        atom_count_path = "/topology/atom_count"
        if not self.contains("topology.spgt.h5", atom_count_path):
            self._validation_problem(
                "topology.spgt.h5 is missing /topology/atom_count"
            )
            return
        atom_count = int(
            self.read_scalar("topology.spgt.h5", atom_count_path)
        )
        for dataset_path in (
            "/atoms/mass",
            "/atoms/charge",
            "/atoms/residue_index",
        ):
            if not self.contains("topology.spgt.h5", dataset_path):
                self._validation_problem(
                    f"topology.spgt.h5 is missing {dataset_path}"
                )
                continue
            shape = self._handles["topology.spgt.h5"][dataset_path].shape
            if shape != (atom_count,):
                raise BundleValidationError(
                    f"topology.spgt.h5 {dataset_path} has shape {shape}, "
                    f"expected ({atom_count},)"
                )

        restart_position = "/particles/all/position/value"
        if self.contains("restart.spgr.h5", restart_position):
            shape = self._handles["restart.spgr.h5"][restart_position].shape
            if len(shape) != 3 or shape[1:] != (atom_count, 3):
                raise BundleValidationError(
                    f"restart.spgr.h5 {restart_position} has shape {shape}, "
                    f"expected (*, {atom_count}, 3)"
                )

    def _validate_finalized_restart(self) -> None:
        status_path = "/parameters/sponge/output/status"
        if not self.has_bundle_file("restart.spgr.h5"):
            return
        if not self.contains("restart.spgr.h5", status_path):
            self._validation_problem(
                f"restart.spgr.h5 is missing {status_path}"
            )
            return
        status = self.read_text("restart.spgr.h5", status_path)
        if status != "finalized":
            self._validation_problem(
                f"restart.spgr.h5 status is {status!r}, expected 'finalized'"
            )

    def _validate_restart_load_policy(self) -> None:
        if not self.has_bundle_file("restart.spgr.h5"):
            return
        policy = (self.case.restart_load or "structural").strip().lower()
        if policy not in {"structural", "dynamic", "protocol", "full"}:
            self._validation_problem(
                "input_h5_restart_load must be structural, dynamic, protocol, "
                f"or full; got {policy!r}"
            )
            return
        has_structural = self.contains(
            "restart.spgr.h5", "/particles/all/position/value"
        )
        if policy == "structural" and not has_structural:
            self._validation_problem(
                "input_h5_restart_load='structural' has no position state"
            )

    def _schema_problem(self, message: str) -> None:
        if self.strict:
            raise BundleSchemaError(message)
        self.warnings.append(message)

    def _validation_problem(self, message: str) -> None:
        if self.strict:
            raise BundleValidationError(message)
        self.warnings.append(message)

    def _require_handle(self, bundle_file: str):
        handle = self._handles.get(bundle_file)
        if handle is None:
            raise BundleValidationError(
                f"bundle artifact is not configured: {bundle_file}"
            )
        return handle


def _decode_h5_text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8")
    if isinstance(value, np.ndarray) and value.shape == ():
        return _decode_h5_text(value.item())
    return str(value)
