"""Legacy surface helpers that patch core runtime types."""

import weakref

import numpy as np

from ._core import (
    Molecule,
    ResidueType,
    add_molecule,
    configure_residue_template_head,
    configure_residue_template_tail,
    get_template_molecule,
    has_template,
    molecule_from_residuetype,
    register_residue_templates_from_mol2_text,
)

_legacy_template_metadata = {}
_legacy_residue_links_override = weakref.WeakKeyDictionary()


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

    @property
    def atoms(self):
        return get_template_molecule(self._name).residues[0].atoms

    def deepcopy(self, name):
        name = str(name)
        _register_template_variant(self._name, name)
        if self.head:
            _legacy_template_metadata.setdefault(name, {})["head"] = self.head
            configure_residue_template_head(name, str(self.head))
        if self.tail:
            _legacy_template_metadata.setdefault(name, {})["tail"] = self.tail
            configure_residue_template_tail(name, str(self.tail))
        return _LegacyResidueTypeHandle(name)

    def omit_atoms(self, atom_names, charge=None):
        del charge  # Legacy API accepts this argument; current migration ignores explicit charge rebalance.
        omit_names = {str(atom_name) for atom_name in atom_names}
        kept_names = [atom.name for atom in self.atoms if atom.name not in omit_names]
        _register_template_variant(self._name, self._name, keep_names=kept_names)
        return self

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
_core_molecule_residue_links = Molecule.residue_links


def _legacy_add_residue_link(self, atom1, atom2):
    pair = [_coerce_atom_index(atom1), _coerce_atom_index(atom2)]
    override = _legacy_residue_links_override.get(self)
    if override is not None:
        if pair not in override:
            override.append(pair)
        return None
    return _core_molecule_add_residue_link(self, pair[0], pair[1])


def _single_residue_mol2_text(template_name, residue_name, keep_names):
    template = get_template_molecule(template_name)
    residue = template.residues[0]
    keep_names = set(keep_names)
    ordered_atoms = [atom for atom in residue.atoms if atom.name in keep_names]
    serial_by_index = {int(atom.index): serial for serial, atom in enumerate(ordered_atoms, start=1)}
    lines = [
        "@<TRIPOS>MOLECULE",
        residue_name,
        f"{len(ordered_atoms)} {sum(1 for atom1, atom2 in template.explicit_bonds if int(atom1) in serial_by_index and int(atom2) in serial_by_index)} 1",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    for serial, atom in enumerate(ordered_atoms, start=1):
        lines.append(
            f"{serial} {atom.name} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f} "
            f"{atom.type} 1 {residue_name} {atom.charge:.6f}"
        )
    lines.append("@<TRIPOS>BOND")
    bond_index = 1
    for atom1, atom2 in template.explicit_bonds:
        atom1 = int(atom1)
        atom2 = int(atom2)
        if atom1 not in serial_by_index or atom2 not in serial_by_index:
            continue
        lines.append(f"{bond_index} {serial_by_index[atom1]} {serial_by_index[atom2]} 1")
        bond_index += 1
    return "\n".join(lines) + "\n"


def _register_template_variant(template_name, residue_name, keep_names=None):
    template = get_template_molecule(template_name)
    residue = template.residues[0]
    if keep_names is None:
        keep_names = [atom.name for atom in residue.atoms]
    text = _single_residue_mol2_text(template_name, residue_name, keep_names)
    register_residue_templates_from_mol2_text(text)
    return residue_name


def _legacy_make_residue_like(value, directly_copy=True):
    del directly_copy  # Core copy policy is currently handled by deepcopying the one-residue molecule.
    if isinstance(value, Molecule):
        if value.residue_count != 1:
            raise TypeError("Residue(...) compatibility only accepts one-residue molecules")
        return value.deepcopy()
    if isinstance(value, ResidueType):
        return molecule_from_residuetype(value)
    if hasattr(value, "name") and has_template(value.name):
        return get_template_molecule(value.name).deepcopy()
    raise TypeError("Residue(...) compatibility expects a one-residue template-like object")


def _legacy_add_residue(self, residue_like):
    add_molecule(self, _legacy_make_residue_like(residue_like))
    return self


def _normalize_residue_link(link):
    if hasattr(link, "atom1") and hasattr(link, "atom2"):
        return [_coerce_atom_index(link.atom1), _coerce_atom_index(link.atom2)]
    if isinstance(link, (list, tuple)) and len(link) == 2:
        return [_coerce_atom_index(link[0]), _coerce_atom_index(link[1])]
    raise TypeError("residue link entries should be (atom1, atom2) pairs or link objects")


def _legacy_get_residue_links(self):
    override = _legacy_residue_links_override.get(self)
    if override is not None:
        return [list(link) for link in override]
    return [list(link) for link in _core_molecule_residue_links.__get__(self, type(self))]


def _legacy_clear_residue_links(self):
    _legacy_residue_links_override[self] = []
    return self


def _legacy_set_residue_links(self, links):
    normalized = []
    seen = set()
    for link in links:
        atom1, atom2 = _normalize_residue_link(link)
        pair = (atom1, atom2)
        if pair in seen:
            continue
        seen.add(pair)
        normalized.append([atom1, atom2])
    _legacy_residue_links_override[self] = normalized
    return self


def _legacy_add_residue_links(self, links):
    current = _legacy_get_residue_links(self)
    current.extend(_normalize_residue_link(link) for link in links)
    return _legacy_set_residue_links(self, current)


def _legacy_get_residue_links_copy(self, copy=True):
    links = _legacy_get_residue_links(self)
    if copy:
        return [list(link) for link in links]
    return links
