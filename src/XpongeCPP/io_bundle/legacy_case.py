"""Discovery and mdin preservation for direct/legacy SPONGE cases."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from .case import parse_mdin_text


_KEY_VALUE_RE = re.compile(
    r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*(?:#.*)?$"
)
_SECTION_RE = re.compile(r"^\s*\[([A-Za-z_][A-Za-z0-9_.]*)\]\s*(?:#.*)?$")


@dataclass(frozen=True)
class LegacyCase:
    """A scanned direct/legacy SPONGE input directory."""

    root: Path
    mdin_path: Path
    mdin_text: str
    commands: dict[str, str]

    @property
    def mode(self) -> str:
        return self.commands.get("mode", "normal").strip().lower()

    def resolve_value_path(self, value: str) -> Path:
        path = Path(value)
        return path if path.is_absolute() else self.root / path

    def resolve_legacy_input_path(self, key: str) -> Path | None:
        value = self.commands.get(key)
        if value:
            return self.resolve_value_path(value)
        if not key.endswith("_in_file"):
            return None
        prefix = self.commands.get("default_in_file_prefix")
        if not prefix:
            return None
        stem = key.removesuffix("_in_file").rstrip("_")
        candidate = self.resolve_value_path(f"{prefix}_{stem}.txt")
        return candidate if candidate.is_file() else None


def render_mdin_without_keys(
    text: str, omit_keys: set[str], append_lines: list[str]
) -> str:
    """Remove normalized legacy bindings while preserving other TOML text."""

    output: list[str] = []
    pending_section: str | None = None
    pending_lines: list[str] = []
    section_has_payload = False
    section: str | None = None

    def flush_section() -> None:
        nonlocal pending_section, pending_lines, section_has_payload
        if pending_section is not None and section_has_payload:
            output.extend(pending_lines)
        pending_section = None
        pending_lines = []
        section_has_payload = False

    for line in text.splitlines():
        section_match = _SECTION_RE.match(line)
        if section_match:
            flush_section()
            section = section_match.group(1).replace(".", "_")
            pending_section = section
            pending_lines = [line]
            section_has_payload = False
            continue
        match = _KEY_VALUE_RE.match(line)
        normalized_key = None
        if match:
            key = match.group(1)
            normalized_key = f"{section}_{key}" if section else key
        if normalized_key in omit_keys:
            continue
        target = pending_lines if pending_section is not None else output
        target.append(line)
        if pending_section is not None and (
            normalized_key is not None or line.strip()
        ):
            section_has_payload = True
    flush_section()
    return "\n".join(output + append_lines).rstrip() + "\n"


def scan_legacy_case(
    case_root: str | Path, mdin: str | Path = "mdin.spg.toml"
) -> LegacyCase:
    """Read a direct/legacy mdin and resolve its case root."""

    root = Path(case_root).resolve()
    mdin_path = Path(mdin)
    if not mdin_path.is_absolute():
        mdin_path = root / mdin_path
    mdin_path = mdin_path.resolve()
    if not mdin_path.is_file():
        raise FileNotFoundError(f"legacy mdin file does not exist: {mdin_path}")
    text = mdin_path.read_text(encoding="utf-8")
    return LegacyCase(root, mdin_path, text, parse_mdin_text(text))
