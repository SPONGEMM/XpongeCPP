"""Discovery and path binding for bundled SPONGE input artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shlex

from .errors import BundleValidationError


_KEY_VALUE_RE = re.compile(
    r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*(?:#.*)?$"
)
_SECTION_RE = re.compile(r"^\s*\[([A-Za-z_][A-Za-z0-9_.]*)\]\s*(?:#.*)?$")
_BUNDLE_BINDINGS = {
    "topology_path": "input_h5_topology_path",
    "protocol_path": "input_h5_protocol_path",
    "restart_path": "input_h5_restart_path",
    "trajectory_path": "input_h5_trajectory_path",
}


def parse_mdin_text(text: str) -> dict[str, str]:
    """Parse the flat key/value subset used by SPONGE input TOML."""

    commands: dict[str, str] = {}
    section: str | None = None
    for line in text.splitlines():
        section_match = _SECTION_RE.match(line)
        if section_match:
            section = section_match.group(1).replace(".", "_")
            continue
        match = _KEY_VALUE_RE.match(line)
        if not match:
            continue
        key, raw_value = match.groups()
        if section:
            key = f"{section}_{key}"
        raw_value = raw_value.strip()
        if not raw_value:
            commands[key] = ""
            continue
        try:
            parts = shlex.split(raw_value, comments=False, posix=True)
        except ValueError:
            parts = []
        commands[key] = (
            parts[0]
            if len(parts) == 1
            else raw_value.strip('"').strip("'")
        )
    return commands


@dataclass(frozen=True)
class BundleCase:
    """Resolved bundled mdin and its HDF5 artifact bindings."""

    root: Path
    mdin_path: Path | None
    mdin_text: str
    commands: dict[str, str]
    topology_path: Path | None
    protocol_path: Path | None
    restart_path: Path | None
    trajectory_path: Path | None = None
    restart_load: str | None = None
    particle_stream: str = "all"

    @property
    def mode(self) -> str:
        return self.commands.get("mode", "normal").strip().lower()

    def path_for_bundle_file(self, bundle_file: str) -> Path | None:
        return {
            "topology.spgt.h5": self.topology_path,
            "protocol.spgp.h5": self.protocol_path,
            "restart.spgr.h5": self.restart_path,
            "trajectory.spg.h5md": self.trajectory_path,
        }.get(bundle_file)


def bundle_case_from_prefix(
    bundle_root: str | Path,
    prefix: str,
    *,
    strict: bool = True,
) -> BundleCase:
    """Bind the canonical three-file input bundle emitted by XpongeCPP."""

    root = Path(bundle_root).resolve()
    paths = {
        "topology_path": root / f"{prefix}_topology.spgt.h5",
        "protocol_path": root / f"{prefix}_protocol.spgp.h5",
        "restart_path": root / f"{prefix}_restart.spgr.h5",
    }
    if strict:
        for field_name, path in paths.items():
            if not path.is_file():
                raise BundleValidationError(
                    f"{field_name.removesuffix('_path')} bundle does not exist: {path}"
                )
    commands = {
        "input_h5_topology_path": str(paths["topology_path"]),
        "input_h5_protocol_path": str(paths["protocol_path"]),
        "input_h5_restart_path": str(paths["restart_path"]),
        "input_h5_restart_load": "structural",
    }
    return BundleCase(
        root=root,
        mdin_path=None,
        mdin_text="",
        commands=commands,
        topology_path=paths["topology_path"],
        protocol_path=paths["protocol_path"],
        restart_path=paths["restart_path"],
        restart_load="structural",
    )


def scan_bundle_case(
    bundle_root: str | Path,
    mdin: str | Path = "mdin.bundled.spg.toml",
    *,
    strict: bool = True,
) -> BundleCase:
    """Read a bundled mdin and resolve its configured HDF5 artifacts."""

    root = Path(bundle_root).resolve()
    mdin_path = Path(mdin)
    if not mdin_path.is_absolute():
        mdin_path = root / mdin_path
    mdin_path = mdin_path.resolve()
    if not mdin_path.is_file():
        raise FileNotFoundError(f"bundled mdin file does not exist: {mdin_path}")

    mdin_text = mdin_path.read_text(encoding="utf-8")
    commands = parse_mdin_text(mdin_text)
    resolved: dict[str, Path | None] = {}
    for field_name, key in _BUNDLE_BINDINGS.items():
        raw_path = commands.get(key)
        if not raw_path:
            resolved[field_name] = None
            continue
        path = Path(raw_path)
        if not path.is_absolute():
            path = mdin_path.parent / path
        path = path.resolve()
        if strict and not path.is_file():
            raise BundleValidationError(f"{key} does not exist: {path}")
        resolved[field_name] = path

    if strict and resolved["topology_path"] is None:
        raise BundleValidationError(
            "bundled mdin does not define input_h5_topology_path"
        )
    if strict and resolved["protocol_path"] is None:
        raise BundleValidationError(
            "bundled mdin does not define input_h5_protocol_path"
        )
    if (
        strict
        and resolved["restart_path"] is None
        and resolved["trajectory_path"] is None
    ):
        raise BundleValidationError(
            "bundled mdin requires input_h5_restart_path or "
            "input_h5_trajectory_path"
        )

    return BundleCase(
        root=root,
        mdin_path=mdin_path,
        mdin_text=mdin_text,
        commands=commands,
        topology_path=resolved["topology_path"],
        protocol_path=resolved["protocol_path"],
        restart_path=resolved["restart_path"],
        trajectory_path=resolved["trajectory_path"],
        restart_load=commands.get("input_h5_restart_load"),
        particle_stream=commands.get(
            "input_h5_trajectory_particle_stream", "all"
        ),
    )
