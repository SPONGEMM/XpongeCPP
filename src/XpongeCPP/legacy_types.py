"""Legacy surface helpers that patch core runtime types."""

import numpy as np

from ._core import Molecule, ResidueType, configure_residue_template_head, configure_residue_template_tail, has_template

_legacy_template_metadata = {}


class _LegacyResidueTypeHandle:
    def __init__(self, name):
        if not has_template(name):
            raise KeyError(f"ResidueType {name!r} is not registered")
        self._name = str(name)

    @property
    def name(self):
        return self._name

    @property
    def head(self):
        return _legacy_template_metadata.get(self._name, {}).get("head")

    @head.setter
    def head(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["head"] = value
        if value:
            configure_residue_template_head(self._name, str(value))

    @property
    def tail(self):
        return _legacy_template_metadata.get(self._name, {}).get("tail")

    @tail.setter
    def tail(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["tail"] = value
        if value:
            configure_residue_template_tail(self._name, str(value))

    def __repr__(self):
        return f"<LegacyResidueTypeHandle {self._name}>"


def _legacy_get_residuetype(name):
    if not has_template(name):
        raise KeyError(f"ResidueType {name!r} is not registered")
    return _LegacyResidueTypeHandle(name)


def _coerce_atom_index(atom):
    if isinstance(atom, (int, np.integer)):
        return int(atom)
    if hasattr(atom, "index"):
        return int(atom.index)
    raise TypeError("atom references should be atom objects or integer indices")


class _AtomIndexProxy:
    def __getitem__(self, atom):
        return _coerce_atom_index(atom)

    def get(self, atom, default=None):
        try:
            return _coerce_atom_index(atom)
        except Exception:
            return default


_core_molecule_add_residue_link = Molecule.add_residue_link


def _legacy_add_residue_link(self, atom1, atom2):
    return _core_molecule_add_residue_link(self, _coerce_atom_index(atom1), _coerce_atom_index(atom2))
