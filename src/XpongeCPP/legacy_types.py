"""Legacy surface helpers that patch core runtime types."""

import weakref

import numpy as np

from ._compat.residuetype import add_template_like, repeat_template_like
from ._core import (
    Molecule,
    ResidueType,
    add_molecule,
    configure_residue_template_head,
    configure_residue_template_tail,
    get_template_molecule,
    has_template,
    molecule_from_residuetype,
    register_residue_type_template,
    registered_template_names,
    register_residue_templates_from_mol2_text,
)

_legacy_template_metadata = {}
_legacy_residue_links_override = weakref.WeakKeyDictionary()
_legacy_dynamic_residue_types = {}


def _remember_dynamic_residuetype(residue_type):
    """Remember writable, Python-created residue types by their public name."""

    _legacy_dynamic_residue_types[str(residue_type.name)] = residue_type
    return residue_type


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
        configure_residue_template_head(
            self._name,
            "" if value is None else str(value),
        )

    @property
    def tail(self):
        return _legacy_template_metadata.get(self._name, {}).get("tail")

    @tail.setter
    def tail(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["tail"] = value
        configure_residue_template_tail(
            self._name,
            "" if value is None else str(value),
        )

    @property
    def head_next(self):
        return _legacy_template_metadata.get(self._name, {}).get("head_next")

    @head_next.setter
    def head_next(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["head_next"] = value

    @property
    def tail_next(self):
        return _legacy_template_metadata.get(self._name, {}).get("tail_next")

    @tail_next.setter
    def tail_next(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["tail_next"] = value

    @property
    def head_length(self):
        return _legacy_template_metadata.get(self._name, {}).get("head_length")

    @head_length.setter
    def head_length(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["head_length"] = value

    @property
    def tail_length(self):
        return _legacy_template_metadata.get(self._name, {}).get("tail_length")

    @tail_length.setter
    def tail_length(self, value):
        _legacy_template_metadata.setdefault(self._name, {})["tail_length"] = value

    @property
    def head_link_conditions(self):
        return _legacy_template_metadata.setdefault(self._name, {}).setdefault(
            "head_link_conditions", []
        )

    @property
    def tail_link_conditions(self):
        return _legacy_template_metadata.setdefault(self._name, {}).setdefault(
            "tail_link_conditions", []
        )

    @property
    def atoms(self):
        dynamic = _legacy_dynamic_residue_types.get(self._name)
        if dynamic is not None:
            return dynamic.atoms
        return get_template_molecule(self._name).residues[0].atoms

    def name2atom(self, name):
        for atom in self.atoms:
            if atom.name == name:
                return atom
        raise KeyError(name)

    def add_atom(self, name, atom_type, x, y, z, charge=0.0, mass=0.0):
        residue_type = _materialize_dynamic_residuetype(self._name)
        return residue_type.add_atom(
            name,
            atom_type,
            x,
            y,
            z,
            charge=charge,
            mass=mass,
        )

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

    def __add__(self, other):
        return add_template_like(self, other)

    def __radd__(self, other):
        return add_template_like(other, self)

    def __mul__(self, count):
        return repeat_template_like(self, count)

    def __rmul__(self, count):
        return repeat_template_like(self, count)

    def __repr__(self):
        return f"<LegacyResidueTypeHandle {self._name}>"


def _legacy_get_residuetype(name):
    dynamic = _legacy_dynamic_residue_types.get(str(name))
    if dynamic is not None:
        return dynamic
    if not has_template(name):
        raise KeyError(f"ResidueType {name!r} is not registered")
    return _LegacyResidueTypeHandle(name)


def _legacy_get_all_residuetypes():
    names = set(registered_template_names())
    names.update(_legacy_dynamic_residue_types)
    return {name: _legacy_get_residuetype(name) for name in sorted(names)}


def _materialize_dynamic_residuetype(name):
    """Return a writable copy of a registered residue template."""

    name = str(name)
    dynamic = _legacy_dynamic_residue_types.get(name)
    if dynamic is not None:
        return dynamic
    if not has_template(name):
        raise KeyError(f"ResidueType {name!r} is not registered")

    template = get_template_molecule(name)
    if template.residue_count != 1:
        raise TypeError(
            "residue template compatibility expects one-residue templates"
        )
    residue = template.residues[0]
    dynamic = _remember_dynamic_residuetype(ResidueType(name))
    for atom in residue.atoms:
        dynamic.add_atom(
            atom.name,
            atom.type,
            atom.x,
            atom.y,
            atom.z,
            charge=atom.charge,
            mass=atom.mass,
        )
    for atom1, atom2 in template.explicit_bonds:
        dynamic.add_connectivity(
            template.atoms[int(atom1)].name,
            template.atoms[int(atom2)].name,
        )
    return dynamic


def _remember_template_connection(
    residue_name,
    position,
    anchor,
    next_atom,
    length,
    conditions=None,
):
    metadata = _legacy_template_metadata.setdefault(str(residue_name), {})
    metadata[str(position)] = anchor
    metadata[f"{position}_next"] = next_atom
    metadata[f"{position}_length"] = length
    if conditions is not None:
        metadata[f"{position}_link_conditions"] = list(conditions)


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
_core_molecule_clear_residue_links = Molecule.clear_residue_links
_core_molecule_residue_links = Molecule.residue_links


def _legacy_add_residue_link(self, atom1, atom2):
    pair = [_coerce_atom_index(atom1), _coerce_atom_index(atom2)]
    _legacy_residue_links_override.pop(self, None)
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


def _dynamic_residuetype_mol2_text(residue_type):
    molecule = molecule_from_residuetype(residue_type)
    residue = molecule.residues[0]
    serial_by_index = {
        int(atom.index): serial
        for serial, atom in enumerate(residue.atoms, start=1)
    }
    lines = [
        "@<TRIPOS>MOLECULE",
        residue_type.name,
        f"{len(residue.atoms)} {len(molecule.explicit_bonds)} 1",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    for serial, atom in enumerate(residue.atoms, start=1):
        lines.append(
            f"{serial} {atom.name} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f} "
            f"{atom.type} 1 {residue_type.name} {atom.charge:.12f}"
        )
    lines.append("@<TRIPOS>BOND")
    for bond_index, (atom1, atom2) in enumerate(
        molecule.explicit_bonds,
        start=1,
    ):
        lines.append(
            f"{bond_index} {serial_by_index[int(atom1)]} "
            f"{serial_by_index[int(atom2)]} 1"
        )
    return "\n".join(lines) + "\n"


def _ensure_dynamic_residuetype_template(name):
    residue_type = _legacy_dynamic_residue_types.get(str(name))
    if residue_type is None:
        return False
    register_residue_type_template(residue_type)
    metadata = _legacy_template_metadata.get(str(name), {})
    if metadata.get("head"):
        configure_residue_template_head(str(name), str(metadata["head"]))
    if metadata.get("tail"):
        configure_residue_template_tail(str(name), str(metadata["tail"]))
    return True


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
    _core_molecule_clear_residue_links(self)
    _legacy_residue_links_override.pop(self, None)
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
    _core_molecule_clear_residue_links(self)
    for atom1, atom2 in normalized:
        _core_molecule_add_residue_link(self, atom1, atom2)
    _legacy_residue_links_override.pop(self, None)
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
