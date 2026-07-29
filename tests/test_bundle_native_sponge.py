"""Direct native-bundle execution gate for a bundle-capable SPONGE runtime."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess

import pytest

from XpongeCPP.io_bundle import convert_legacy_to_bundle

from test_bundle_converter import _legacy_case, _read_mdout_first_frame


def _run(executable: Path, case_dir: Path, mdin: str):
    completed = subprocess.run(
        [str(executable), "-mdin", mdin],
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert completed.returncode == 0, (
        f"bundle-capable SPONGE failed in {case_dir}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    return _read_mdout_first_frame(case_dir / "mdout.txt")


def test_native_bundle_matches_raw_sponge_zeroth_frame(tmp_path):
    configured = os.environ.get("SPONGE_BUNDLE_EXECUTABLE")
    if not configured:
        pytest.skip(
            "set SPONGE_BUNDLE_EXECUTABLE to a runtime supporting input_h5_*"
        )
    executable = Path(configured).expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        pytest.fail(
            "SPONGE_BUNDLE_EXECUTABLE is not executable: "
            f"{executable}"
        )

    raw_dir = _legacy_case(tmp_path)
    converted_root = tmp_path / "converted"
    convert_legacy_to_bundle(raw_dir, converted_root)
    bundle_dir = converted_root / "bundle"

    raw = _run(executable, raw_dir, "mdin.spg.toml")
    bundled = _run(executable, bundle_dir, "mdin.bundled.spg.toml")

    assert bundled.keys() == raw.keys()
    for key in raw:
        assert bundled[key] == pytest.approx(
            raw[key], abs=1.0e-5, rel=1.0e-7
        ), key
