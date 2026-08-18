"""Shared parmchk2 input helpers for GAFF and GAFF2."""

from __future__ import annotations

import os
from importlib import import_module
from tempfile import TemporaryDirectory

from ... import (
    Molecule,
    Save_Mol2,
    get_template_molecule,
    has_template,
    molecule_from_residuetype,
)


def import_xpongelib():
    try:
        xlib = import_module("XpongeLib")
    except ImportError as exc:
        raise ImportError(
            "parmchk2 requires the external XpongeLib runtime "
            "(install mokda-xpongelib / XpongeLib first)"
        ) from exc
    if not hasattr(xlib, "_parmchk2"):
        raise ImportError("Installed XpongeLib does not expose _parmchk2")
    return xlib


def coerce_parmchk2_input(ifname):
    if isinstance(ifname, (str, os.PathLike)):
        return str(ifname), None
    tempdir = TemporaryDirectory()
    tempfile = os.path.join(tempdir.name, "temp.mol2")
    if isinstance(ifname, Molecule):
        Save_Mol2(ifname, tempfile)
        return tempfile, tempdir
    if (
        hasattr(ifname, "atom_count")
        and hasattr(ifname, "bond_count")
        and hasattr(ifname, "name")
    ):
        Save_Mol2(molecule_from_residuetype(ifname), tempfile)
        return tempfile, tempdir
    if hasattr(ifname, "name") and has_template(ifname.name):
        Save_Mol2(get_template_molecule(ifname.name), tempfile)
        return tempfile, tempdir
    tempdir.cleanup()
    raise TypeError(
        "parmchk2 expects a mol2 path or a template-like molecule object"
    )


__all__ = [
    "coerce_parmchk2_input",
    "import_xpongelib",
]
