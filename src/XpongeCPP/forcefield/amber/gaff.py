"""Amber GAFF compatibility helpers."""

import os
from importlib import import_module
from tempfile import TemporaryDirectory

from ... import (
    Molecule,
    Save_Mol2,
    get_template_molecule,
    has_template,
    molecule_from_residuetype,
    register_amber_parmdat_file,
)
from . import data_path
from . import load_parameters_from_frcmod

register_amber_parmdat_file(str(data_path("gaff.dat")))


def _import_xpongelib():
    try:
        xlib = import_module("XpongeLib")
    except ImportError as exc:
        raise ImportError(
            "gaff.parmchk2_gaff requires the external XpongeLib runtime "
            "(install mokda-xpongelib / XpongeLib first)"
        ) from exc
    if not hasattr(xlib, "_parmchk2"):
        raise ImportError("Installed XpongeLib does not expose _parmchk2")
    return xlib


def _coerce_parmchk2_input(ifname):
    if isinstance(ifname, (str, os.PathLike)):
        return str(ifname), None
    tempdir = TemporaryDirectory()
    tempfile = os.path.join(tempdir.name, "temp.mol2")
    if isinstance(ifname, Molecule):
        Save_Mol2(ifname, tempfile)
        return tempfile, tempdir
    if hasattr(ifname, "atom_count") and hasattr(ifname, "bond_count") and hasattr(ifname, "name"):
        Save_Mol2(molecule_from_residuetype(ifname), tempfile)
        return tempfile, tempdir
    if hasattr(ifname, "name") and has_template(ifname.name):
        Save_Mol2(get_template_molecule(ifname.name), tempfile)
        return tempfile, tempdir
    tempdir.cleanup()
    raise TypeError("parmchk2_gaff expects a mol2 path or a template-like molecule object")


def parmchk2_gaff(ifname, ofname, direct_load=True, keep=True):
    """Generate frcmod parameters with legacy Xponge-compatible semantics."""
    xlib = _import_xpongelib()
    mol2_path, tempdir = _coerce_parmchk2_input(ifname)
    try:
        datapath = os.path.dirname(xlib.__file__)
        xlib._parmchk2(mol2_path, "mol2", str(ofname), datapath, 0, 1, 1)
        if direct_load:
            load_parameters_from_frcmod(ofname, prefix=False)
        if not keep:
            os.remove(ofname)
    finally:
        if tempdir is not None:
            tempdir.cleanup()
