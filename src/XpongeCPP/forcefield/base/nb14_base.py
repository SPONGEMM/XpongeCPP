"""Legacy nb14_base compatibility helpers."""

from __future__ import annotations


class _NB14TypeRecord:
    def __init__(self, name, kLJ, kee):
        self.name = name
        self.kLJ = float(kLJ)
        self.kee = float(kee)


class _NB14ForceEntity:
    def __init__(self, atoms, nb14_type, name=None):
        self.atoms = list(atoms)
        self.type = nb14_type
        self.name = nb14_type.name if name is None else name
        self.kLJ = nb14_type.kLJ
        self.kee = nb14_type.kee


class NB14Type:
    """Minimal executable compatibility for legacy 1-4 nonbonded types."""

    topology_matrix = [[1, -4], [1, 1]]
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
    def entity(cls, atoms, nb14_type, name=None):
        return _NB14ForceEntity(atoms, nb14_type, name=name)

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
            if len(values) < 3:
                continue
            entry = _NB14TypeRecord(values[0], values[1], values[2])
            cls._types[cls._norm(values[0])] = entry
        return cls


__all__ = ["NB14Type"]
