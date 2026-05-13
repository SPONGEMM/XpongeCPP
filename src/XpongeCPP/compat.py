"""Optional compatibility helpers for legacy Xponge calling style."""

from __future__ import annotations

import inspect

from . import (
    Molecule,
    get_template_molecule,
    registered_template_names,
    save_gro,
    save_mol2,
    save_pdb,
    save_sponge_input,
)


def _save_pdb_method(self, filename, write_cryst1=True):
    return save_pdb(self, filename, write_cryst1)


def _save_mol2_method(self, filename):
    return save_mol2(self, filename)


def _save_gro_method(self, filename):
    return save_gro(self, filename)


def _save_sponge_input_method(self, prefix=None, dirname="."):
    return save_sponge_input(self, "" if prefix is None else prefix, dirname)


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


def install_template_globals(namespace=None, template_names=None, overwrite=False):
    """Inject registered template molecules into a Python namespace."""
    if namespace is None:
        frame = inspect.currentframe()
        assert frame is not None and frame.f_back is not None
        namespace = frame.f_back.f_globals
    names = registered_template_names() if template_names is None else list(template_names)
    for name in names:
        if not overwrite and name in namespace:
            continue
        namespace[name] = get_template_molecule(name)
    return namespace


def enable_legacy_namespace(namespace=None, template_names=None, overwrite=False):
    """Enable common legacy Xponge conveniences in one call."""
    install_molecule_io_methods()
    return install_template_globals(namespace=namespace, template_names=template_names, overwrite=overwrite)


Enable_Legacy_Namespace = enable_legacy_namespace
Install_Template_Globals = install_template_globals
Install_Molecule_IO_Methods = install_molecule_io_methods

install_molecule_io_methods()

__all__ = [
    "enable_legacy_namespace",
    "Enable_Legacy_Namespace",
    "install_template_globals",
    "Install_Template_Globals",
    "install_molecule_io_methods",
    "Install_Molecule_IO_Methods",
]
