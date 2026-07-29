import os
import random
import sys
from copy import deepcopy
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA_DIR = REPO_ROOT / "tests" / "data"
DATA_1KV2_DIR = TEST_DATA_DIR / "1kv2"


def pytest_collection_modifyitems(config, items):
    """Provide dependency-free, reproducible order randomization for CI."""
    del config
    seed_text = os.environ.get("XPONGECPP_TEST_ORDER_SEED")
    if seed_text is None:
        return
    try:
        seed = int(seed_text)
    except ValueError as exc:
        raise pytest.UsageError(
            "XPONGECPP_TEST_ORDER_SEED must be an integer"
        ) from exc
    random.Random(seed).shuffle(items)


@pytest.fixture(autouse=True)
def isolate_forcefield_registries():
    """Restore all process-global force-field state after every test."""
    import XpongeCPP
    from XpongeCPP.forcefield.amber import _forcefield_family
    from XpongeCPP.helper import AtomType

    core_snapshot = XpongeCPP._core._snapshot_forcefield_registries()
    active_forcefields = dict(_forcefield_family._ACTIVE_FORCEFIELDS)
    atom_types = deepcopy(AtomType._types)
    atom_parameters = dict(AtomType._parameters)
    atom_property_units = dict(AtomType._property_units)
    module_names = set(sys.modules)
    try:
        yield
    finally:
        XpongeCPP._core._restore_forcefield_registries(core_snapshot)
        _forcefield_family._ACTIVE_FORCEFIELDS.clear()
        _forcefield_family._ACTIVE_FORCEFIELDS.update(active_forcefields)
        AtomType._types.clear()
        AtomType._types.update(atom_types)
        AtomType._parameters.clear()
        AtomType._parameters.update(atom_parameters)
        AtomType._property_units.clear()
        AtomType._property_units.update(atom_property_units)
        prefixes = (
            "XpongeCPP.forcefield.",
            "XpongeCPP.data.",
            "Xponge.forcefield.",
            "Xponge.data.",
        )
        for name in sorted(set(sys.modules) - module_names, reverse=True):
            if name.startswith(prefixes):
                module = sys.modules.pop(name, None)
                parent_name, _, child_name = name.rpartition(".")
                parent = sys.modules.get(parent_name)
                if (
                    module is not None
                    and parent is not None
                    and getattr(parent, child_name, None) is module
                ):
                    delattr(parent, child_name)


def original_xponge_repo() -> Path:
    configured = os.environ.get("XPONGE_REFERENCE_REPO")
    if configured:
        return Path(configured)
    candidates = [
        REPO_ROOT.parent / "Xponge-origin",
        REPO_ROOT.parent / "Xponge",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def optional_1kv2_baseline_dir() -> Path | None:
    configured = os.environ.get("XPONGECPP_1KV2_BASELINE_DIR")
    if not configured:
        candidates = [
            TEST_DATA_DIR / "1kv2_xponge_baseline",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None
    return Path(configured)
