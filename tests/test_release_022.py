"""Release-contract checks for the Xponge-origin 1.7b9 alignment."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
import tomllib

import XpongeCPP


ROOT = Path(__file__).resolve().parents[1]


def test_release_version_and_compatibility_target():
    metadata = tomllib.loads(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    assert metadata["version"] == XpongeCPP.__version__ == "0.2.2"
    assert "Xponge-origin 1.7b9" in metadata["description"]
    assert metadata["requires-python"] == ">=3.10,<3.13"


def test_cli_version_matches_imported_version():
    completed = subprocess.run(
        [sys.executable, "-m", "XpongeCPP", "-v"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.strip() == XpongeCPP.__version__
    assert completed.stderr == ""


def test_metal_assignment_is_the_only_metal_workflow_namespace():
    removed_name = "MC" + "PB"
    assert XpongeCPP.metal_assignment.__name__ == "XpongeCPP.metal_assignment"
    assert not hasattr(XpongeCPP, removed_name)
    assert importlib.util.find_spec(
        "XpongeCPP." + removed_name.lower()
    ) is None


def test_plain_import_is_silent():
    completed = subprocess.run(
        [sys.executable, "-c", "import XpongeCPP"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout == ""
    assert completed.stderr == ""
