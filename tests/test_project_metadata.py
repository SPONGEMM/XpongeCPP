import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_runtime_dependencies_track_xponge_origin_distribution_set():
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert metadata["project"]["dependencies"] == [
        "numpy<2",
        "h5py>=3.11",
        "geometric>=1.1",
        "matplotlib>=3.10.8",
        "MDAnalysis>=2.9.0",
        "PubChemPy>=1.0.5",
        "rdkit>=2025.9.3",
        "pyscf>=2.11.0; platform_system != 'Windows'",
        "mokda-xpongelib>=1.2.5.0; platform_system == 'Darwin' or platform_machine == 'x86_64' or platform_machine == 'AMD64'",
    ]


def test_distribution_version_matches_runtime_version():
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    source = (ROOT / "src" / "XpongeCPP" / "__init__.py").read_text(encoding="utf-8")
    assert f'__version__ = "{metadata["project"]["version"]}"' in source
