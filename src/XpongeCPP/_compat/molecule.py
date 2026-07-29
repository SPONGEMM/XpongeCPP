"""Legacy molecule-surface helpers centralized under the compat package."""

import weakref

from .. import Molecule, save_gro, save_pdb, save_sponge_input


_BUILT_STATE = weakref.WeakKeyDictionary()


def _save_pdb_method(self, filename, write_cryst1=True):
    return save_pdb(self, filename, write_cryst1)


def _save_mol2_method(self, filename=None):
    from ..process import Save_Mol2

    return Save_Mol2(self, filename)


def _save_gro_method(self, filename):
    return save_gro(self, filename)


def _save_sponge_input_method(self, prefix=None, dirname="."):
    from ..process import Save_SPONGE_Input

    return Save_SPONGE_Input(self, prefix, dirname)


def _get_atoms_method(self):
    """Origin-compatible no-op cache hook for native eager atom storage."""

    return self.atoms


def _get_built(self):
    return _BUILT_STATE.get(self, True)


def _set_built(self, value):
    _BUILT_STATE[self] = bool(value)


def install_molecule_io_methods():
    """Attach instance-style save methods to ``Molecule`` for legacy scripts."""
    method_map = {
        "save_pdb": _save_pdb_method,
        "Save_PDB": _save_pdb_method,
        "save_mol2": _save_mol2_method,
        "Save_Mol2": _save_mol2_method,
        "save_gro": _save_gro_method,
        "Save_GRO": _save_gro_method,
        "save_sponge_input": _save_sponge_input_method,
        "Save_SPONGE_Input": _save_sponge_input_method,
        "get_atoms": _get_atoms_method,
        "Get_Atoms": _get_atoms_method,
    }
    for name, func in method_map.items():
        if not hasattr(Molecule, name):
            setattr(Molecule, name, func)
    if not hasattr(Molecule, "_all"):
        Molecule._all = {}
    if not hasattr(Molecule, "built"):
        Molecule.built = property(_get_built, _set_built)
