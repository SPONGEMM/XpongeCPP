"""Convert validated SPONGE bundles into direct/legacy input files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import numpy as np

from .case import BundleCase, bundle_case_from_prefix, scan_bundle_case
from .errors import (
    BundleCapabilityError,
    BundleExportError,
    BundleValidationError,
)
from .legacy_materializer import LegacyMaterializer, LegacyPayload
from .manifest import ManifestEntry, ReverseConversionManifest
from .reader import BundleReader


_TOPOLOGY = "topology.spgt.h5"
_RESTART = "restart.spgr.h5"
_UNSUPPORTED_ROOTS = (
    "/forcefield/virtual_atom",
    "/forcefield/urey_bradley",
    "/forcefield/bond_soft",
    "/forcefield/custom_force",
    "/forcefield/gb",
    "/forcefield/subsys_division",
    "/forcefield/cmap",
    "/forcefield/lj_soft_core",
)
_H5_INPUT_KEYS = {
    "input_h5_topology_path",
    "input_h5_protocol_path",
    "input_h5_restart_path",
    "input_h5_restart_load",
    "input_h5_trajectory_path",
    "input_h5_trajectory_particle_stream",
    "default_in_file_prefix",
}


def _text_values(values) -> list[str]:
    result = []
    for value in np.asarray(values).reshape(-1):
        if isinstance(value, (bytes, np.bytes_)):
            result.append(bytes(value).decode("utf-8"))
        else:
            result.append(str(value))
    return result


def _fixed_rows(rows, precision: int = 6) -> str:
    return "".join(
        " ".join(
            str(int(value))
            if isinstance(value, (int, np.integer))
            else f"{float(value):.{precision}f}"
            for value in row
        )
        + "\n"
        for row in rows
    )


def _box_dimensions(edges) -> np.ndarray:
    cell = np.asarray(edges, dtype=np.float64).reshape(3, 3)
    lengths = np.linalg.norm(cell, axis=1)
    if np.any(lengths <= 0):
        raise BundleExportError("bundle box edges must have positive lengths")

    def angle(lhs, rhs, lhs_length, rhs_length):
        cosine = np.dot(lhs, rhs) / (lhs_length * rhs_length)
        return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))

    return np.asarray(
        [
            *lengths,
            angle(cell[1], cell[2], lengths[1], lengths[2]),
            angle(cell[0], cell[2], lengths[0], lengths[2]),
            angle(cell[0], cell[1], lengths[0], lengths[1]),
        ]
    )


class BundleToLegacyConverter:
    """Materialize the lossless core of a native XpongeCPP input bundle."""

    def __init__(
        self,
        case: BundleCase,
        output_dir: str | Path,
        *,
        prefix: str | None = None,
        strict: bool = True,
        overwrite: bool = False,
    ):
        self.case = case
        self.output_dir = Path(output_dir).resolve()
        self.prefix = prefix or "input"
        self.strict = strict
        self.materializer = LegacyMaterializer(
            self.output_dir, overwrite=overwrite
        )
        self.manifest = ReverseConversionManifest(
            bundle_root=str(case.root),
            output_root=str(self.output_dir),
            mode=case.mode,
        )
        self._bindings: dict[str, str] = {}

    def convert(self, *, dry_run: bool = False) -> ReverseConversionManifest:
        """Plan, validate, and optionally write direct SPONGE inputs."""

        with BundleReader(self.case, strict=self.strict) as reader:
            self.manifest.warnings.extend(reader.warnings)
            self._validate_capabilities(reader)
            exporters: tuple[tuple[str, Callable[[BundleReader], str]], ...] = (
                ("residue", self._export_residue),
                ("resname", self._export_resname),
                ("atom_name", self._export_atom_name),
                ("atom_type_name", self._export_atom_type_name),
                ("mass", self._export_mass),
                ("charge", self._export_charge),
                ("coordinate", self._export_coordinate),
                ("LJ", self._export_lj),
                ("bond", self._export_bond),
                ("angle", self._export_angle),
                ("dihedral", self._export_dihedral),
                ("exclude", self._export_exclude),
                ("nb14", self._export_nb14),
            )
            if reader.contains(_TOPOLOGY, "/forcefield/improper/atoms"):
                exporters += (("improper_dihedral", self._export_improper),)
            for key, exporter in exporters:
                self._plan(key, exporter(reader), self._source_for(key))
            self._restore_sidecars(reader)

        self._plan_mdin()
        if not dry_run:
            manifest_path = self.output_dir / "manifest.bundle_to_legacy.json"
            manifest_payload = json.dumps(
                self.manifest.to_dict(), indent=2, sort_keys=True
            ) + "\n"
            self.materializer.plan(
                LegacyPayload("__manifest__", manifest_payload),
                manifest_path.name,
            )
            self.materializer.write_all()
        else:
            self.materializer.validate_targets()
        return self.manifest

    def _problem(self, message: str) -> None:
        if self.strict:
            raise BundleCapabilityError(message)
        self.manifest.warnings.append(message)

    def _validate_capabilities(self, reader: BundleReader) -> None:
        for root in _UNSUPPORTED_ROOTS:
            if reader.contains(_TOPOLOGY, root):
                if root == "/forcefield/lj_soft_core":
                    raise BundleCapabilityError(
                        "legacy conversion does not yet support typed data at "
                        f"{root}"
                    )
                self._problem(
                    f"legacy conversion does not yet support typed data at {root}"
                )
        params = np.asarray(reader.read(_TOPOLOGY, "/forcefield/nb14/params"))
        if params.ndim != 2 or params.shape[1] != 2:
            raise BundleCapabilityError(
                "legacy nb14 conversion requires params shape (N, 2); "
                f"got {params.shape}"
            )
        if reader.has_bundle_file("protocol.spgp.h5"):
            for key in ("cv_count", "restraint_count"):
                path = f"/protocol/{key}"
                if reader.contains("protocol.spgp.h5", path):
                    value = int(
                        reader.read_scalar("protocol.spgp.h5", path)
                    )
                    if value:
                        self._problem(
                            f"legacy conversion does not yet support {key}={value}"
                        )
        try:
            charge_unit = reader.read_attribute(
                _TOPOLOGY, "/atoms/charge", "unit"
            )
        except BundleValidationError:
            self._problem("/atoms/charge must declare unit='Amber'")
        else:
            if charge_unit != "Amber":
                self._problem(
                    f"/atoms/charge unit is {charge_unit!r}, expected 'Amber'"
                )

    def _source_for(self, key: str) -> str:
        if key == "coordinate":
            return str(self.case.restart_path)
        return str(self.case.topology_path)

    def _plan(
        self, key: str, payload: str | bytes, source_path: str
    ) -> Path:
        filename = f"{self.prefix}_{key}.txt"
        target = self.materializer.plan(LegacyPayload(key, payload), filename)
        self._bindings[f"{key}_in_file"] = filename
        self.manifest.add(
            ManifestEntry(
                key=key,
                source_path=source_path,
                target_path=str(target),
            )
        )
        return target

    def _restore_sidecars(self, reader: BundleReader) -> None:
        for bundle_file in (_TOPOLOGY, "protocol.spgp.h5", _RESTART):
            if not reader.has_bundle_file(bundle_file):
                continue
            for key, source in reader.read_legacy_sidecars(bundle_file).items():
                if f"{key}_in_file" in self._bindings:
                    continue
                target = self._plan(
                    key,
                    source.read_bytes(),
                    str(source),
                )
                self.manifest.entries[-1] = ManifestEntry(
                    key=key,
                    source_path=str(source),
                    target_path=str(target),
                    status="sidecar_restored",
                )

    def _plan_mdin(self) -> None:
        lines = []
        for key, value in self.case.commands.items():
            if key not in _H5_INPUT_KEYS and not key.endswith("_in_file"):
                lines.append(f'{key} = "{value}"')
        for key, value in sorted(self._bindings.items()):
            lines.append(f'{key} = "{value}"')
        payload = "\n".join(lines) + "\n"
        target = self.materializer.plan(
            LegacyPayload("run_mdin", payload), "mdin.legacy.spg.toml"
        )
        self.manifest.generated_mdin = str(target)
        self.manifest.add(
            ManifestEntry(
                key="run_mdin",
                source_path=str(self.case.mdin_path or self.case.root),
                target_path=str(target),
                status="mdin_binding_generated",
            )
        )

    def _read_required(self, reader: BundleReader, path: str):
        if not reader.contains(_TOPOLOGY, path):
            raise BundleValidationError(f"{_TOPOLOGY} is missing {path}")
        return reader.read(_TOPOLOGY, path)

    def _export_residue(self, reader: BundleReader) -> str:
        atom_count = int(reader.read_scalar(_TOPOLOGY, "/topology/atom_count"))
        offsets = np.asarray(
            self._read_required(reader, "/residues/atom_offset"), dtype=np.int64
        )
        if offsets.ndim != 1 or len(offsets) < 1 or offsets[0] != 0:
            raise BundleExportError("invalid /residues/atom_offset")
        counts = np.diff(offsets)
        if offsets[-1] != atom_count or np.any(counts < 0):
            raise BundleExportError(
                "/residues/atom_offset does not cover all atoms"
            )
        return (
            f"{atom_count} {len(counts)}\n"
            + "".join(f"{int(value)}\n" for value in counts)
        )

    def _export_text_vector(
        self, reader: BundleReader, path: str, expected: int
    ) -> str:
        values = _text_values(self._read_required(reader, path))
        if len(values) != expected:
            raise BundleExportError(
                f"{path} has length {len(values)}, expected {expected}"
            )
        return f"{expected}\n" + "".join(f"{value}\n" for value in values)

    def _export_resname(self, reader: BundleReader) -> str:
        offsets = self._read_required(reader, "/residues/atom_offset")
        return self._export_text_vector(
            reader, "/parameters/xponge/residues/name", len(offsets) - 1
        )

    def _export_atom_name(self, reader: BundleReader) -> str:
        count = int(reader.read_scalar(_TOPOLOGY, "/topology/atom_count"))
        return self._export_text_vector(
            reader, "/parameters/xponge/atoms/name", count
        )

    def _export_atom_type_name(self, reader: BundleReader) -> str:
        count = int(reader.read_scalar(_TOPOLOGY, "/topology/atom_count"))
        return self._export_text_vector(
            reader, "/parameters/xponge/atoms/type_name", count
        )

    def _export_mass(self, reader: BundleReader) -> str:
        values = np.asarray(self._read_required(reader, "/atoms/mass")).reshape(-1)
        return f"{len(values)}\n" + "".join(
            f"{float(value):.3f}\n" for value in values
        )

    def _export_charge(self, reader: BundleReader) -> str:
        values = np.asarray(
            self._read_required(reader, "/atoms/charge")
        ).reshape(-1)
        return f"{len(values)}\n" + "".join(
            f"{float(value):.6f}\n" for value in values
        )

    def _export_coordinate(self, reader: BundleReader) -> str:
        positions = np.asarray(
            reader.read(_RESTART, "/particles/all/position/value")
        )
        edges = np.asarray(
            reader.read(_RESTART, "/particles/all/box/edges/value")
        )
        if positions.ndim != 3 or positions.shape[-1] != 3:
            raise BundleExportError(
                f"restart positions have invalid shape {positions.shape}"
            )
        if edges.ndim != 3 or edges.shape[1:] != (3, 3):
            raise BundleExportError(
                f"restart box edges have invalid shape {edges.shape}"
            )
        frame = positions[-1]
        box = _box_dimensions(edges[-1])
        return (
            f"{len(frame)}\n"
            + _fixed_rows(frame)
            + " ".join(f"{value:.6f}" for value in box)
            + "\n"
        )

    def _export_lj(self, reader: BundleReader) -> str:
        type_count = int(
            reader.read_scalar(_TOPOLOGY, "/forcefield/lj/atom_type_count")
        )
        pair_a = np.asarray(
            self._read_required(reader, "/forcefield/lj/pair_A_12")
        ).reshape(-1)
        pair_b = np.asarray(
            self._read_required(reader, "/forcefield/lj/pair_B_6")
        ).reshape(-1)
        atom_types = np.asarray(
            self._read_required(reader, "/forcefield/lj/type"), dtype=np.int64
        ).reshape(-1)
        pair_count = type_count * (type_count + 1) // 2
        if len(pair_a) != pair_count or len(pair_b) != pair_count:
            raise BundleExportError("LJ triangular table has invalid length")
        rows_a, rows_b = [], []
        cursor = 0
        for width in range(1, type_count + 1):
            rows_a.append(pair_a[cursor : cursor + width])
            rows_b.append(pair_b[cursor : cursor + width])
            cursor += width
        render = lambda rows: "".join(  # noqa: E731
            " ".join(f"{float(value):.6e}" for value in row) + " \n"
            for row in rows
        )
        return (
            f"{len(atom_types)} {type_count}\n\n"
            + render(rows_a)
            + "\n"
            + render(rows_b)
            + "\n"
            + "".join(f"{int(value)}\n" for value in atom_types)
        )

    def _export_interaction(
        self,
        reader: BundleReader,
        root: str,
        scalar_paths: tuple[str, ...],
    ) -> str:
        atoms = np.asarray(self._read_required(reader, f"{root}/atoms"))
        scalars = [
            np.asarray(self._read_required(reader, f"{root}/{name}")).reshape(-1)
            for name in scalar_paths
        ]
        if atoms.ndim != 2 or any(len(values) != len(atoms) for values in scalars):
            raise BundleExportError(f"{root} interaction arrays have mismatched shapes")
        rows = [
            tuple(int(value) for value in atom_row)
            + tuple(
                int(values[index])
                if name == "periodicity"
                else float(values[index])
                for name, values in zip(scalar_paths, scalars)
            )
            for index, atom_row in enumerate(atoms)
        ]
        rows.sort()
        return f"{len(rows)}\n" + _fixed_rows(rows)

    def _export_bond(self, reader: BundleReader) -> str:
        return self._export_interaction(
            reader, "/forcefield/bond", ("k", "r0")
        )

    def _export_angle(self, reader: BundleReader) -> str:
        return self._export_interaction(
            reader, "/forcefield/angle", ("k", "theta0")
        )

    def _export_dihedral(self, reader: BundleReader) -> str:
        return self._export_interaction(
            reader,
            "/forcefield/dihedral",
            ("periodicity", "k", "phi0"),
        )

    def _export_improper(self, reader: BundleReader) -> str:
        return self._export_interaction(
            reader, "/forcefield/improper", ("pk", "phi0")
        )

    def _export_exclude(self, reader: BundleReader) -> str:
        offsets = np.asarray(
            self._read_required(reader, "/topology/exclusions/offset"),
            dtype=np.int64,
        )
        values = np.asarray(
            self._read_required(reader, "/topology/exclusions/list"),
            dtype=np.int64,
        ).reshape(-1)
        atom_count = int(reader.read_scalar(_TOPOLOGY, "/topology/atom_count"))
        if offsets.shape != (atom_count + 1,) or offsets[0] != 0:
            raise BundleExportError("exclusion offsets have invalid shape")
        if offsets[-1] != len(values) or np.any(np.diff(offsets) < 0):
            raise BundleExportError("exclusion offsets do not cover the list")
        lines = [f"{atom_count} {len(values)}"]
        for index in range(atom_count):
            row = values[offsets[index] : offsets[index + 1]]
            lines.append(
                f"{len(row)} " + " ".join(str(int(value)) for value in row)
            )
        return "\n".join(lines) + "\n"

    def _export_nb14(self, reader: BundleReader) -> str:
        atoms = np.asarray(
            self._read_required(reader, "/forcefield/nb14/atoms")
        )
        params = np.asarray(
            self._read_required(reader, "/forcefield/nb14/params")
        )
        if atoms.ndim != 2 or atoms.shape[1] != 2:
            raise BundleExportError(
                f"/forcefield/nb14/atoms has invalid shape {atoms.shape}"
            )
        if params.shape != (len(atoms), 2):
            raise BundleExportError(
                f"/forcefield/nb14/params has invalid shape {params.shape}"
            )
        rows = [
            (
                int(atom_row[0]),
                int(atom_row[1]),
                float(param_row[0]),
                float(param_row[1]),
            )
            for atom_row, param_row in zip(atoms, params)
        ]
        return f"{len(rows)}\n" + _fixed_rows(rows)


def convert_bundle_to_legacy(
    bundle_root: str | Path,
    output_dir: str | Path,
    *,
    mdin: str | Path | None = "mdin.bundled.spg.toml",
    prefix: str | None = None,
    strict: bool = True,
    overwrite: bool = False,
    dry_run: bool = False,
) -> ReverseConversionManifest:
    """Convert a scanned bundle or a native ``prefix`` bundle to legacy files."""

    root = Path(bundle_root).resolve()
    mdin_path = root / mdin if mdin is not None else None
    if mdin_path is not None and mdin_path.is_file():
        case = scan_bundle_case(root, mdin_path, strict=strict)
    elif prefix is not None:
        case = bundle_case_from_prefix(root, prefix, strict=strict)
    else:
        raise FileNotFoundError(
            "bundled mdin was not found; pass prefix=... for a native "
            "XpongeCPP bundle"
        )
    return BundleToLegacyConverter(
        case,
        output_dir,
        prefix=prefix,
        strict=strict,
        overwrite=overwrite,
    ).convert(dry_run=dry_run)


__all__ = ["BundleToLegacyConverter", "convert_bundle_to_legacy"]
