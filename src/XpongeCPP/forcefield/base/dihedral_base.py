"""Legacy dihedral_base compatibility helpers."""

from __future__ import annotations

from itertools import permutations


class _DihedralTypeRecord:
    def __init__(self, name, k, phi0, periodicity):
        self.name = name
        self.k = float(k)
        self.phi0 = float(phi0)
        self.periodicity = int(periodicity)


class _DihedralForceEntity:
    def __init__(self, atoms, dihedral_type, name=None):
        self.atoms = list(atoms)
        self.type = dihedral_type
        self.name = dihedral_type.name if name is None else name
        self.k = dihedral_type.k
        self.phi0 = dihedral_type.phi0
        self.periodicity = dihedral_type.periodicity
        # Old Xponge proper torsions can behave like a multi-term container.
        self.ks = [self.k]
        self.phi0s = [self.phi0]
        self.periodicitys = [self.periodicity]
        self.multiple_numbers = 1


class _BaseDihedralType:
    _types = {}

    @staticmethod
    def _norm(name):
        return str(name).strip().upper()

    @classmethod
    def get_type(cls, name):
        return cls._types[cls._norm(name)]

    @classmethod
    def getType(cls, name):
        return cls.get_type(name)

    @classmethod
    def Get_Type(cls, name):
        return cls.get_type(name)

    @classmethod
    def get_all_types(cls):
        return dict(cls._types)

    @classmethod
    def Get_All_Types(cls):
        return cls.get_all_types()

    @classmethod
    def Clear_Type(cls):
        cls._types.clear()

    @classmethod
    def entity(cls, atoms, dihedral_type, name=None):
        return _DihedralForceEntity(atoms, dihedral_type, name=name)

    @classmethod
    def _store(cls, name, k, phi0, periodicity):
        entry = _DihedralTypeRecord(name, k, phi0, periodicity)
        cls._types[cls._norm(name)] = entry
        return entry

    @classmethod
    def New_From_String(cls, text):
        lines = [line.strip() for line in str(text).splitlines() if line.strip()]
        if not lines:
            return cls
        headers = lines[0].split()
        if not headers or headers[0].lower() != "name":
            return cls
        for line in lines[1:]:
            values = line.split()
            if len(values) < 4:
                continue
            cls._store(values[0], values[1], values[2], values[3])
        return cls


class ProperType(_BaseDihedralType):
    """Minimal executable compatibility for legacy proper dihedral types."""

    _types = {}


class ImproperType(_BaseDihedralType):
    """Minimal executable compatibility for legacy improper dihedral types."""

    _types = {}
    topology_matrix = [[1, 3, 2, 3], [1, 1, 2, 3], [1, 1, 1, 2], [1, 1, 1, 1]]

    @staticmethod
    def same_force(atom_list):
        """Match the first-wave old Xponge permutation semantics."""
        temp = []
        if isinstance(atom_list, str):
            atom_list_temp = [atom.strip() for atom in atom_list.split("-")]
            center_atom = atom_list_temp.pop(2)
            for atom_permutation in permutations(atom_list_temp):
                atom_permutation = list(atom_permutation)
                atom_permutation.insert(2, center_atom)
                temp.append("-".join(atom_permutation))
        else:
            atom_list_temp = list(atom_list)
            center_atom = atom_list_temp.pop(2)
            for atom_permutation in permutations(atom_list_temp):
                atom_permutation = list(atom_permutation)
                atom_permutation.insert(2, center_atom)
                temp.append(atom_permutation)
        return temp


__all__ = ["ProperType", "ImproperType"]
