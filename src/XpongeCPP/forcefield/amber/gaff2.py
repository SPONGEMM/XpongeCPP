"""Register Amber GAFF2 parameters from packaged data."""

from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("small_molecule", "gaff2")

import os

from ... import AtomType, implemented_gaff2_assign_types, register_amber_parmdat_file
from . import data_path
from . import load_parameters_from_frcmod
from ._parmchk2 import (
    coerce_parmchk2_input,
    filter_mixed_gaff_mol2,
    import_xpongelib,
)

register_amber_parmdat_file(str(data_path("gaff2.dat")))
AtomType.New_From_String("\n".join(implemented_gaff2_assign_types()))


def parmchk2_gaff2(ifname, ofname, direct_load=True, keep=True):
    """Generate GAFF2 frcmod parameters with Xponge-compatible semantics."""
    xlib = import_xpongelib()
    mol2_path, tempdir = coerce_parmchk2_input(ifname)
    filtered_tempdir = None
    try:
        mol2_path, filtered_tempdir = filter_mixed_gaff_mol2(mol2_path)
        datapath = os.path.dirname(xlib.__file__)
        xlib._parmchk2(mol2_path, "mol2", str(ofname), datapath, 0, 1, 2)
        if direct_load:
            load_parameters_from_frcmod(ofname, prefix=False)
        if not keep:
            os.remove(ofname)
    finally:
        if filtered_tempdir is not None:
            filtered_tempdir.cleanup()
        if tempdir is not None:
            tempdir.cleanup()
