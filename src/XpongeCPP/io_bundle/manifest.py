"""Machine-readable manifests for bundle conversion."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class ManifestEntry:
    """One materialized legacy payload."""

    key: str
    source_path: str
    target_path: str
    status: str = "typed_exported"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ReverseConversionManifest:
    """Summary returned by bundle-to-legacy conversion."""

    schema: str = "xponge.bundle_to_legacy.manifest"
    schema_version: int = 1
    bundle_root: str = ""
    output_root: str = ""
    mode: str = "normal"
    entries: list[ManifestEntry] = field(default_factory=list)
    generated_mdin: str | None = None
    warnings: list[str] = field(default_factory=list)

    def add(self, entry: ManifestEntry) -> None:
        self.entries.append(entry)

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema": self.schema,
            "schema_version": self.schema_version,
            "bundle_root": self.bundle_root,
            "output_root": self.output_root,
            "mode": self.mode,
            "entries": [entry.to_dict() for entry in self.entries],
        }
        if self.generated_mdin is not None:
            data["generated_mdin"] = self.generated_mdin
        if self.warnings:
            data["warnings"] = list(self.warnings)
        return data
