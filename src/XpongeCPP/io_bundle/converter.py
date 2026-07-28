"""Convert ordinary direct/legacy SPONGE inputs into v2 bundles."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import tempfile
from uuid import uuid4

import numpy as np

from .bundle_builder import (
    canonical_dataset_hash,
    write_dataset,
    write_string,
)
from .errors import (
    BundleCapabilityError,
    BundleConflictError,
    BundleValidationError,
)
from .legacy_case import (
    LegacyCase,
    render_mdin_without_keys,
    scan_legacy_case,
)
from .manifest import ConversionManifest, ManifestEntry


_REQUIRED_KEYS = (
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
_OPTIONAL_KEYS = ("improper_dihedral",)
_INPUT_KEYS = {
    f"{key}_in_file" for key in (*_REQUIRED_KEYS, *_OPTIONAL_KEYS)
}
_BUNDLE_BINDINGS = (
    'input_h5_topology_path = "topology.spgt.h5"',
    'input_h5_protocol_path = "protocol.spgp.h5"',
    'input_h5_restart_path = "restart.spgr.h5"',
    'input_h5_restart_load = "structural"',
)


def _tokens(path: Path) -> list[str]:
    try:
        return path.read_text(encoding="utf-8").split()
    except UnicodeDecodeError as exc:
        raise BundleValidationError(
            f"legacy input must be UTF-8 text: {path}"
        ) from exc


def _integer(token: str, path: Path, label: str) -> int:
    try:
        value = int(token)
    except ValueError as exc:
        raise BundleValidationError(
            f"{path} has invalid integer for {label}: {token!r}"
        ) from exc
    return value


def _floating(token: str, path: Path, label: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise BundleValidationError(
            f"{path} has invalid float for {label}: {token!r}"
        ) from exc
    if not np.isfinite(value):
        raise BundleValidationError(
            f"{path} has non-finite value for {label}: {token!r}"
        )
    return value


def _counted_vector(
    path: Path, *, dtype, label: str
) -> np.ndarray:
    tokens = _tokens(path)
    if not tokens:
        raise BundleValidationError(f"{path} is empty")
    count = _integer(tokens[0], path, f"{label} count")
    if count < 0 or len(tokens) != count + 1:
        raise BundleValidationError(
            f"{path} declares {count} {label} values but has "
            f"{len(tokens) - 1}"
        )
    parser = _integer if np.issubdtype(np.dtype(dtype), np.integer) else _floating
    return np.asarray(
        [
            parser(token, path, f"{label}[{index}]")
            for index, token in enumerate(tokens[1:])
        ],
        dtype=dtype,
    )


def _counted_text(path: Path, label: str) -> np.ndarray:
    tokens = _tokens(path)
    if not tokens:
        raise BundleValidationError(f"{path} is empty")
    count = _integer(tokens[0], path, f"{label} count")
    values = tokens[1:]
    if count < 0 or len(values) != count:
        raise BundleValidationError(
            f"{path} declares {count} {label} values but has {len(values)}"
        )
    return np.asarray(values, dtype=object)


def _interaction(
    path: Path,
    *,
    atom_columns: int,
    parameter_columns: tuple[tuple[str, np.dtype], ...],
) -> dict[str, np.ndarray]:
    tokens = _tokens(path)
    if not tokens:
        raise BundleValidationError(f"{path} is empty")
    count = _integer(tokens[0], path, "interaction count")
    width = atom_columns + len(parameter_columns)
    values = tokens[1:]
    if count < 0 or len(values) != count * width:
        raise BundleValidationError(
            f"{path} declares {count} rows of width {width}, "
            f"got {len(values)} values"
        )
    rows = [values[index : index + width] for index in range(0, len(values), width)]
    atoms = np.asarray(
        [
            [
                _integer(value, path, f"atom column {column}")
                for column, value in enumerate(row[:atom_columns])
            ]
            for row in rows
        ],
        dtype=np.int32,
    ).reshape(count, atom_columns)
    result = {"atoms": atoms}
    for offset, (name, dtype) in enumerate(parameter_columns):
        parser = _integer if np.issubdtype(dtype, np.integer) else _floating
        result[name] = np.asarray(
            [
                parser(row[atom_columns + offset], path, name)
                for row in rows
            ],
            dtype=dtype,
        )
    result["count"] = np.asarray(count, dtype=np.int64)
    return result


def _box_edges(dimensions: np.ndarray) -> np.ndarray:
    a, b, c, alpha_deg, beta_deg, gamma_deg = dimensions
    alpha, beta, gamma = np.radians([alpha_deg, beta_deg, gamma_deg])
    sin_gamma = np.sin(gamma)
    if abs(sin_gamma) < 1e-7:
        raise BundleValidationError("legacy box gamma produces a singular cell")
    edge_y_x = b * np.cos(gamma)
    edge_y_y = b * sin_gamma
    edge_z_x = c * np.cos(beta)
    edge_z_y = c * (
        np.cos(alpha) - np.cos(beta) * np.cos(gamma)
    ) / sin_gamma
    z_squared = c * c - edge_z_x * edge_z_x - edge_z_y * edge_z_y
    if z_squared < -1e-5:
        raise BundleValidationError("legacy box angles produce an invalid cell")
    return np.asarray(
        [
            [a, 0.0, 0.0],
            [edge_y_x, edge_y_y, 0.0],
            [edge_z_x, edge_z_y, np.sqrt(max(0.0, z_squared))],
        ],
        dtype=np.float32,
    )


class LegacyToBundleConverter:
    """Convert the ordinary XpongeCPP raw-output subset into a v2 bundle."""

    def __init__(self, case: LegacyCase, output_dir: str | Path):
        self.case = case
        self.output_dir = Path(output_dir).resolve()
        self.bundle_dir = self.output_dir / "bundle"
        self.manifest = ConversionManifest(
            case_root=str(case.root), mode=case.mode
        )
        self._sources: dict[str, Path] = {}

    def convert(self, *, dry_run: bool = False) -> ConversionManifest:
        """Validate the whole legacy case, then publish its bundle."""

        if self.case.mode not in {"normal", "minimization", "md"}:
            raise BundleCapabilityError(
                f"legacy-to-bundle conversion does not support mode "
                f"{self.case.mode!r}"
            )
        self._resolve_sources()
        topology, atom_count = self._parse_topology()
        restart = self._parse_restart(atom_count)
        self._record_manifest()
        final_mdin = self.bundle_dir / "mdin.bundled.spg.toml"
        self.manifest.bundled_mdin = str(final_mdin)
        if dry_run:
            self._validate_targets()
            return self.manifest

        self._validate_targets()
        self.output_dir.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(
                prefix=f".{self.output_dir.name}.",
                dir=self.output_dir.parent,
            )
        )
        try:
            stage_bundle = staging / "bundle"
            stage_bundle.mkdir()
            self._write_bundle(stage_bundle, topology, restart, atom_count)
            mdin = render_mdin_without_keys(
                self.case.mdin_text,
                _INPUT_KEYS | {"default_in_file_prefix"},
                list(_BUNDLE_BINDINGS),
            )
            (stage_bundle / "mdin.bundled.spg.toml").write_text(
                mdin, encoding="utf-8"
            )
            (staging / "manifest.json").write_text(
                json.dumps(
                    self.manifest.to_dict(), indent=2, sort_keys=True
                )
                + "\n",
                encoding="utf-8",
            )
            os.replace(staging, self.output_dir)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        return self.manifest

    def _validate_targets(self) -> None:
        if self.output_dir.exists():
            raise BundleConflictError(
                f"bundle output directory already exists: {self.output_dir}"
            )

    def _resolve_sources(self) -> None:
        for key in _REQUIRED_KEYS:
            command_key = f"{key}_in_file"
            path = self.case.resolve_legacy_input_path(command_key)
            if path is None or not path.is_file():
                raise BundleValidationError(
                    f"legacy case is missing required {command_key}"
                )
            self._sources[key] = path.resolve()
        for key in _OPTIONAL_KEYS:
            path = self.case.resolve_legacy_input_path(f"{key}_in_file")
            if path is not None:
                if not path.is_file():
                    raise BundleValidationError(
                        f"legacy input does not exist: {path}"
                    )
                self._sources[key] = path.resolve()

        unsupported = []
        for command_key in self.case.commands:
            if (
                command_key.endswith("_in_file")
                and command_key not in _INPUT_KEYS
                and self.case.resolve_legacy_input_path(command_key) is not None
            ):
                unsupported.append(command_key)
        if unsupported:
            raise BundleCapabilityError(
                "legacy-to-bundle conversion does not yet support: "
                + ", ".join(sorted(unsupported))
            )

    def _parse_topology(self) -> tuple[dict[str, np.ndarray], int]:
        datasets: dict[str, np.ndarray] = {}
        residue_tokens = _tokens(self._sources["residue"])
        if len(residue_tokens) < 2:
            raise BundleValidationError("residue input is missing its header")
        atom_count = _integer(
            residue_tokens[0], self._sources["residue"], "atom count"
        )
        residue_count = _integer(
            residue_tokens[1], self._sources["residue"], "residue count"
        )
        counts = np.asarray(
            [
                _integer(value, self._sources["residue"], "residue atom count")
                for value in residue_tokens[2:]
            ],
            dtype=np.int64,
        )
        if counts.shape != (residue_count,) or np.any(counts < 0):
            raise BundleValidationError(
                "residue input count does not match its residue rows"
            )
        if int(counts.sum()) != atom_count:
            raise BundleValidationError(
                "residue atom counts do not sum to the declared atom count"
            )
        offsets = np.concatenate(
            [np.asarray([0], dtype=np.int64), np.cumsum(counts)]
        )
        residue_index = np.repeat(
            np.arange(residue_count, dtype=np.int32), counts
        )
        datasets["/residues/atom_offset"] = offsets
        datasets["/atoms/residue_index"] = residue_index

        mass = _counted_vector(
            self._sources["mass"], dtype=np.float32, label="mass"
        )
        charge = _counted_vector(
            self._sources["charge"], dtype=np.float32, label="charge"
        )
        if len(mass) != atom_count or len(charge) != atom_count:
            raise BundleValidationError(
                "mass and charge counts must match the residue atom count"
            )
        datasets["/atoms/mass"] = mass
        datasets["/atoms/charge"] = charge

        text_specs = (
            ("resname", "/parameters/xponge/residues/name", residue_count),
            ("atom_name", "/parameters/xponge/atoms/name", atom_count),
            (
                "atom_type_name",
                "/parameters/xponge/atoms/type_name",
                atom_count,
            ),
        )
        for key, dataset_path, expected in text_specs:
            values = _counted_text(self._sources[key], key)
            if len(values) != expected:
                raise BundleValidationError(
                    f"{key} count is {len(values)}, expected {expected}"
                )
            datasets[dataset_path] = values

        interactions = (
            (
                "bond",
                "/forcefield/bond",
                2,
                (("k", np.dtype("float32")), ("r0", np.dtype("float32"))),
            ),
            (
                "angle",
                "/forcefield/angle",
                3,
                (
                    ("k", np.dtype("float32")),
                    ("theta0", np.dtype("float32")),
                ),
            ),
            (
                "dihedral",
                "/forcefield/dihedral",
                4,
                (
                    ("periodicity", np.dtype("int32")),
                    ("k", np.dtype("float32")),
                    ("phi0", np.dtype("float32")),
                ),
            ),
        )
        for key, root, atom_columns, parameter_columns in interactions:
            parsed = _interaction(
                self._sources[key],
                atom_columns=atom_columns,
                parameter_columns=parameter_columns,
            )
            for name, values in parsed.items():
                datasets[f"{root}/{name}"] = values

        if "improper_dihedral" in self._sources:
            parsed = _interaction(
                self._sources["improper_dihedral"],
                atom_columns=4,
                parameter_columns=(
                    ("pk", np.dtype("float32")),
                    ("phi0", np.dtype("float32")),
                ),
            )
            for name, values in parsed.items():
                datasets[f"/forcefield/improper/{name}"] = values

        self._parse_nb14(datasets)
        self._parse_exclusions(datasets, atom_count)
        self._parse_lj(datasets, atom_count)
        self._validate_atom_indices(datasets, atom_count)
        return datasets, atom_count

    def _parse_nb14(self, datasets: dict[str, np.ndarray]) -> None:
        parsed = _interaction(
            self._sources["nb14"],
            atom_columns=2,
            parameter_columns=(
                ("k_lj", np.dtype("float32")),
                ("k_ee", np.dtype("float32")),
            ),
        )
        datasets["/forcefield/nb14/atoms"] = parsed["atoms"]
        datasets["/forcefield/nb14/params"] = np.column_stack(
            [parsed["k_lj"], parsed["k_ee"]]
        ).astype(np.float32, copy=False)
        datasets["/forcefield/nb14/count"] = parsed["count"]

    def _parse_exclusions(
        self, datasets: dict[str, np.ndarray], atom_count: int
    ) -> None:
        path = self._sources["exclude"]
        tokens = _tokens(path)
        if len(tokens) < 2:
            raise BundleValidationError("exclude input is missing its header")
        declared_atoms = _integer(tokens[0], path, "exclude atom count")
        declared_values = _integer(tokens[1], path, "exclude value count")
        if declared_atoms != atom_count:
            raise BundleValidationError(
                f"exclude atom count is {declared_atoms}, expected {atom_count}"
            )
        cursor = 2
        offsets = [0]
        values: list[int] = []
        for atom_index in range(atom_count):
            if cursor >= len(tokens):
                raise BundleValidationError(
                    f"exclude input ends before atom {atom_index}"
                )
            row_count = _integer(tokens[cursor], path, "exclude row count")
            cursor += 1
            row = tokens[cursor : cursor + row_count]
            if len(row) != row_count:
                raise BundleValidationError("exclude row is truncated")
            values.extend(
                _integer(value, path, "excluded atom") for value in row
            )
            cursor += row_count
            offsets.append(len(values))
        if cursor != len(tokens) or len(values) != declared_values:
            raise BundleValidationError(
                "exclude payload does not match its declared total"
            )
        datasets["/topology/exclusions/offset"] = np.asarray(
            offsets, dtype=np.int64
        )
        datasets["/topology/exclusions/list"] = np.asarray(
            values, dtype=np.int32
        )

    def _parse_lj(
        self, datasets: dict[str, np.ndarray], atom_count: int
    ) -> None:
        path = self._sources["LJ"]
        tokens = _tokens(path)
        if len(tokens) < 2:
            raise BundleValidationError("LJ input is missing its header")
        declared_atoms = _integer(tokens[0], path, "LJ atom count")
        type_count = _integer(tokens[1], path, "LJ type count")
        pair_count = type_count * (type_count + 1) // 2
        expected = 2 + pair_count * 2 + atom_count
        if declared_atoms != atom_count or len(tokens) != expected:
            raise BundleValidationError(
                "LJ payload dimensions do not match its header"
            )
        cursor = 2
        pair_a = np.asarray(
            [
                _floating(value, path, "LJ A")
                for value in tokens[cursor : cursor + pair_count]
            ],
            dtype=np.float32,
        )
        cursor += pair_count
        pair_b = np.asarray(
            [
                _floating(value, path, "LJ B")
                for value in tokens[cursor : cursor + pair_count]
            ],
            dtype=np.float32,
        )
        cursor += pair_count
        atom_type = np.asarray(
            [
                _integer(value, path, "LJ atom type")
                for value in tokens[cursor:]
            ],
            dtype=np.int32,
        )
        if np.any(atom_type < 0) or np.any(atom_type >= type_count):
            raise BundleValidationError("LJ atom type is outside its table")
        datasets["/forcefield/lj/atom_type_count"] = np.asarray(
            type_count, dtype=np.int32
        )
        datasets["/forcefield/lj/pair_A_12"] = pair_a
        datasets["/forcefield/lj/pair_B_6"] = pair_b
        datasets["/forcefield/lj/type"] = atom_type

    def _validate_atom_indices(
        self, datasets: dict[str, np.ndarray], atom_count: int
    ) -> None:
        for path, values in datasets.items():
            if not path.endswith("/atoms"):
                continue
            if np.any(values < 0) or np.any(values >= atom_count):
                raise BundleValidationError(
                    f"{path} contains an out-of-range atom index"
                )
        exclusions = datasets["/topology/exclusions/list"]
        if np.any(exclusions < 0) or np.any(exclusions >= atom_count):
            raise BundleValidationError(
                "exclusion list contains an out-of-range atom index"
            )

    def _parse_restart(
        self, atom_count: int
    ) -> dict[str, np.ndarray]:
        path = self._sources["coordinate"]
        tokens = _tokens(path)
        expected = 1 + atom_count * 3 + 6
        if len(tokens) != expected:
            raise BundleValidationError(
                f"{path} has {len(tokens)} fields, expected {expected}"
            )
        declared_atoms = _integer(tokens[0], path, "coordinate atom count")
        if declared_atoms != atom_count:
            raise BundleValidationError(
                f"coordinate atom count is {declared_atoms}, "
                f"expected {atom_count}"
            )
        values = np.asarray(
            [
                _floating(token, path, "coordinate")
                for token in tokens[1:]
            ],
            dtype=np.float64,
        )
        positions = values[: atom_count * 3].reshape(1, atom_count, 3)
        edges = _box_edges(values[atom_count * 3 :]).reshape(1, 3, 3)
        return {
            "/particles/all/position/value": positions.astype(np.float32),
            "/particles/all/box/edges/value": edges,
            "/particles/all/step": np.asarray([0], dtype=np.int64),
            "/particles/all/time": np.asarray([0.0], dtype=np.float64),
        }

    def _record_manifest(self) -> None:
        bundle_file_by_key = {
            "coordinate": "restart.spgr.h5",
        }
        for key, source in sorted(self._sources.items()):
            bundle_file = bundle_file_by_key.get(key, "topology.spgt.h5")
            self.manifest.add(
                ManifestEntry(
                    key=key,
                    source_path=str(source),
                    target_path=str(self.bundle_dir / bundle_file),
                    status="typed_converted",
                    bundle_file=bundle_file,
                )
            )

    def _write_bundle(
        self,
        root: Path,
        topology: dict[str, np.ndarray],
        restart: dict[str, np.ndarray],
        atom_count: int,
    ) -> None:
        import h5py

        identity_uuid = str(uuid4())
        topology_hash = canonical_dataset_hash(
            "topology.spgt.h5", topology
        )
        atom_order_hash = canonical_dataset_hash(
            "topology.spgt.h5",
            topology,
            path_prefixes=("/atoms/", "/residues/"),
        )
        forcefield_hash = canonical_dataset_hash(
            "topology.spgt.h5",
            topology,
            path_prefixes=("/forcefield/", "/manybody/", "/qc/"),
        )
        protocol_hash = canonical_dataset_hash("protocol.spgp.h5", {})
        state_hash = canonical_dataset_hash(
            "restart.spgr.h5",
            restart,
            path_prefixes=("/particles/", "/parameters/restart/"),
        )

        topology_path = root / "topology.spgt.h5"
        with h5py.File(topology_path, "w") as handle:
            for path, value in topology.items():
                write_dataset(
                    handle,
                    path,
                    value,
                    string=np.asarray(value).dtype.kind in {"O", "U", "S"},
                )
            handle["/atoms/charge"].attrs["unit"] = "Amber"
            self._write_schema(
                handle, "sponge.topology.h5", identity_uuid
            )
            write_string(handle, "/topology/atom_order_hash", atom_order_hash)
            write_string(handle, "/topology/topology_hash", topology_hash)
            write_string(handle, "/topology/forcefield_hash", forcefield_hash)
            write_dataset(
                handle,
                "/topology/atom_count",
                np.asarray(atom_count, dtype=np.int64),
            )

        protocol_path = root / "protocol.spgp.h5"
        with h5py.File(protocol_path, "w") as handle:
            self._write_schema(
                handle, "sponge.protocol.h5", identity_uuid
            )
            write_string(
                handle,
                "/protocol/topology_compatibility/topology_hash",
                topology_hash,
            )
            write_string(handle, "/identity/content_hash", protocol_hash)
            write_dataset(
                handle,
                "/protocol/cv_count",
                np.asarray(0, dtype=np.int64),
            )
            write_dataset(
                handle,
                "/protocol/restraint_count",
                np.asarray(0, dtype=np.int64),
            )

        restart_path = root / "restart.spgr.h5"
        with h5py.File(restart_path, "w") as handle:
            for group in (
                "/h5md",
                "/h5md/creator",
                "/run",
                "/particles/all",
                "/parameters/restart",
                "/parameters/restart/rng_state",
                "/parameters/restart/integrator_state",
                "/parameters/restart/thermostat",
                "/parameters/restart/barostat",
                "/parameters/restart/protocol_sidecars",
                "/parameters/restart/bias",
                "/parameters/restart/bias/sits",
                "/parameters/restart/bias/meta",
            ):
                handle.require_group(group)
            for path, value in restart.items():
                write_dataset(handle, path, value)
            handle["/h5md"].attrs["version"] = np.asarray(
                [1, 1], dtype=np.int32
            )
            handle["/h5md/creator"].attrs["name"] = "XpongeCPP"
            handle["/h5md/creator"].attrs["version"] = "legacy-to-bundle"
            handle["/particles/all/time"].attrs["unit"] = "ps"
            handle["/particles/all/position/value"].attrs["unit"] = "Angstrom"
            handle["/particles/all/box/edges/value"].attrs["unit"] = "Angstrom"
            box = handle["/particles/all/box"]
            box.attrs["dimension"] = np.int32(3)
            box.attrs["boundary"] = np.asarray(
                ["periodic", "periodic", "periodic"],
                dtype=h5py.string_dtype("utf-8"),
            )
            handle["/particles/all/position/step"] = handle[
                "/particles/all/step"
            ]
            handle["/particles/all/position/time"] = handle[
                "/particles/all/time"
            ]
            handle["/particles/all/box/edges/step"] = handle[
                "/particles/all/step"
            ]
            handle["/particles/all/box/edges/time"] = handle[
                "/particles/all/time"
            ]
            self._write_schema(handle, "sponge.restart.h5", identity_uuid)
            write_string(handle, "/run/topology_hash", topology_hash)
            write_string(handle, "/run/atom_order_hash", atom_order_hash)
            write_string(
                handle, "/run/producer_protocol_hash", protocol_hash
            )
            write_string(handle, "/run/state_hash", state_hash)
            write_string(
                handle,
                "/parameters/sponge/output/status",
                "finalized",
            )
            write_dataset(
                handle,
                "/parameters/sponge/output/frame_count",
                np.asarray([1], dtype=np.int64),
            )
            write_dataset(
                handle,
                "/parameters/sponge/output/last_complete_step",
                np.asarray([0], dtype=np.int64),
            )
            write_dataset(
                handle,
                "/parameters/sponge/output/last_complete_time",
                np.asarray([0.0], dtype=np.float64),
            )
            write_dataset(
                handle,
                "/run/current_step",
                np.asarray([0], dtype=np.int64),
            )
            write_dataset(
                handle,
                "/run/current_time",
                np.asarray([0.0], dtype=np.float64),
            )
            write_dataset(
                handle,
                "/parameters/sponge/output/particle_streams",
                np.asarray(["all"], dtype=object),
                string=True,
            )
            write_string(handle, "/run/state_type", "restart")

    @staticmethod
    def _write_schema(handle, schema_name: str, identity_uuid: str) -> None:
        write_string(handle, "/schema/name", schema_name)
        write_string(handle, "/schema/version", "sponge.input.v2")
        write_string(handle, "/parameters/sponge/schema/name", schema_name)
        write_string(
            handle, "/parameters/sponge/schema/version", "sponge.input.v2"
        )
        write_string(handle, "/identity/uuid", identity_uuid)


def convert_legacy_to_bundle(
    case_root: str | Path,
    output_dir: str | Path,
    *,
    mdin: str | Path = "mdin.spg.toml",
    dry_run: bool = False,
) -> ConversionManifest:
    """Convert a direct/legacy SPONGE case into canonical bundle files."""

    case = scan_legacy_case(case_root, mdin)
    return LegacyToBundleConverter(case, output_dir).convert(dry_run=dry_run)


__all__ = ["LegacyToBundleConverter", "convert_legacy_to_bundle"]
