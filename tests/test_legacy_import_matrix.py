import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


CORE_SCOPE_LEGACY_IMPORTS = [
    "Xponge",
    "Xponge.analysis",
    "Xponge.analysis.md_analysis",
    "Xponge.analysis.sasa",
    "Xponge.analysis.wham",
    "Xponge.assign",
    "Xponge.build",
    "Xponge.forcefield.amber.bsc1",
    "Xponge.forcefield.amber.ff14sb",
    "Xponge.forcefield.amber.ff19sb",
    "Xponge.forcefield.amber.gaff",
    "Xponge.forcefield.amber.ol3",
    "Xponge.forcefield.amber.rsff2c",
    "Xponge.forcefield.amber.tip3p",
    "Xponge.forcefield.amber.tip4pew",
    "Xponge.forcefield.base.angle_base",
    "Xponge.forcefield.base.cmap_base",
    "Xponge.forcefield.base.charge_base",
    "Xponge.forcefield.base.dihedral_base",
    "Xponge.forcefield.base.lj_base",
    "Xponge.forcefield.base.mass_base",
    "Xponge.forcefield.base.nb14_base",
    "Xponge.forcefield.base.virtual_atom_base",
    "Xponge.forcefield.charmm.charmm27",
    "Xponge.forcefield.charmm.charmm36",
    "Xponge.forcefield.charmm.tip3p_charmm",
    "Xponge.forcefield.opls.oplsaam",
    "Xponge.forcefield.special.gb",
    "Xponge.forcefield.special.fep",
    "Xponge.forcefield.special.min",
    "Xponge.forcefield.sw.mw",
    "Xponge.helper",
    "Xponge.helper.cv",
    "Xponge.helper.file",
    "Xponge.helper.gromacs",
    "Xponge.helper.math",
    "Xponge.helper.namespace",
    "Xponge.load",
    "Xponge.mdrun",
    "Xponge.process",
    "Xponge.tools",
    "Xponge.tools.unittests",
]


def _run_import(module_name):
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    return subprocess.run(
        [
            sys.executable,
            "-c",
            f"import importlib; importlib.import_module({module_name!r})",
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
    )


def _public_forcefield_imports():
    package_root = ROOT / "src" / "XpongeCPP" / "forcefield"
    modules = set()
    for path in package_root.rglob("*.py"):
        relative = path.relative_to(ROOT / "src").with_suffix("")
        parts = list(relative.parts)
        if parts[-1] == "__init__":
            parts.pop()
        if any(part.startswith("_") for part in parts[2:]):
            continue
        native = ".".join(parts)
        modules.add(native)
        modules.add(native.replace("XpongeCPP", "Xponge", 1))
    return sorted(modules)


@pytest.mark.parametrize("module_name", CORE_SCOPE_LEGACY_IMPORTS)
def test_core_scope_legacy_import_matrix(module_name):
    result = _run_import(module_name)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("module_name", _public_forcefield_imports())
def test_public_forcefield_import_matrix(module_name):
    result = _run_import(module_name)
    assert result.returncode == 0, result.stderr
