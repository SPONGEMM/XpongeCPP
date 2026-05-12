"""Force-field compatibility namespace."""

from importlib import resources
from pathlib import Path


def package_data_path(*parts):
    return resources.files("XpongeCPP").joinpath("data", *parts)


def repository_reference_forcefield_path(*parts):
    repo_root = Path(__file__).resolve().parents[3]
    reference_root = repo_root / "third_party" / "xponge_reference_forcefield"
    if reference_root.exists():
        return reference_root.joinpath(*parts)
    return package_data_path("reference_forcefield", *parts)
