"""Legacy molecule-surface helpers centralized under the compat package."""

from .. import Molecule, save_gro, save_pdb, save_sponge_input


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
    }
    for name, func in method_map.items():
        if not hasattr(Molecule, name):
            setattr(Molecule, name, func)
