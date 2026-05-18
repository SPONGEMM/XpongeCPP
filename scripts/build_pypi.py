#!/usr/bin/env python3
"""Build, validate, and smoke-test a local PyPI distribution."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import venv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist"
BUILD = ROOT / "build"


def run(cmd: list[str], *, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=ROOT, env=env)


def venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def latest_wheel() -> Path:
    wheels = sorted(DIST.glob("*.whl"))
    if not wheels:
        raise RuntimeError("No wheel was produced in dist/")
    return wheels[-1]


def main() -> None:
    shutil.rmtree(DIST, ignore_errors=True)
    shutil.rmtree(BUILD, ignore_errors=True)

    run([sys.executable, "-m", "pip", "install", "--upgrade", "build", "twine"])
    run([sys.executable, "-m", "build"])
    run([sys.executable, "-m", "twine", "check", "dist/*"])

    with tempfile.TemporaryDirectory(prefix="xpongecpp-wheel-smoke-") as tmp:
        venv_dir = Path(tmp) / "venv"
        venv.EnvBuilder(with_pip=True).create(venv_dir)
        py = venv_python(venv_dir)
        run([str(py), "-m", "pip", "install", "--upgrade", "pip"])
        run([str(py), "-m", "pip", "install", str(latest_wheel())])
        smoke = (
            "import XpongeCPP, Xponge; "
            "import Xponge.forcefield.amber.ff19sb; "
            "from Xponge.forcefield.special import gb; "
            "print(XpongeCPP.__version__)"
        )
        run([str(py), "-c", smoke])

    print("Build, twine check, and wheel smoke test all passed.")


if __name__ == "__main__":
    main()
