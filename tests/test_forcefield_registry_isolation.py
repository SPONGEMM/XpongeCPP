import importlib
import sys

import pytest

import XpongeCPP
from XpongeCPP.forcefield.amber import _forcefield_family
from XpongeCPP.helper import AtomType


def _unload_forcefield_module(name):
    sys.modules.pop(name, None)
    _forcefield_family._ACTIVE_FORCEFIELDS.clear()


def test_registry_snapshot_restores_amber_templates_and_python_atom_types():
    snapshot = XpongeCPP._core._snapshot_forcefield_registries()
    template_names = set(XpongeCPP.registered_template_names())
    atom_type_names = set(AtomType.get_all_types())

    importlib.import_module("XpongeCPP.forcefield.amber.ff19sb")
    importlib.import_module("XpongeCPP.forcefield.amber.gaff2")
    assert XpongeCPP.has_template("NHYP")
    assert set(AtomType.get_all_types()) != atom_type_names

    XpongeCPP._core._restore_forcefield_registries(snapshot)
    AtomType._types.clear()
    AtomType._types.update({name: AtomType(name) for name in atom_type_names})
    assert set(XpongeCPP.registered_template_names()) == template_names


@pytest.mark.parametrize(
    ("first", "second", "family", "first_name", "second_name"),
    [
        ("ff14sb", "ff19sb", "protein", "ff14sb", "ff19sb"),
        ("gaff", "gaff2", "small_molecule", "gaff", "gaff2"),
    ],
)
def test_registry_snapshot_allows_in_process_forcefield_transition(
    first, second, family, first_name, second_name
):
    snapshot = XpongeCPP._core._snapshot_forcefield_registries()
    first_module = f"XpongeCPP.forcefield.amber.{first}"
    second_module = f"XpongeCPP.forcefield.amber.{second}"

    importlib.import_module(first_module)
    assert _forcefield_family.get_active_forcefield(family) == first_name

    XpongeCPP._core._restore_forcefield_registries(snapshot)
    _unload_forcefield_module(first_module)
    importlib.import_module(second_module)
    assert _forcefield_family.get_active_forcefield(family) == second_name
