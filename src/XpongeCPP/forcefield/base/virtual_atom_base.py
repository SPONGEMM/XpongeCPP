"""Legacy virtual_atom_base compatibility helpers."""

from __future__ import annotations

from XpongeCPP.gromacs import GlobalSetting


class _VirtualType2Record:
    def __init__(self, name, atom0, atom1, atom2, k1, k2):
        self.name = name
        self.atom0 = int(atom0)
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.k1 = float(k1)
        self.k2 = float(k2)


class VirtualType2:
    """Minimal executable compatibility for legacy type-2 virtual atoms."""

    _types = {}

    @staticmethod
    def _norm(name):
        return str(name).strip().upper()

    @classmethod
    def get_type(cls, name):
        return cls._types[cls._norm(name)]

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
    def New_From_String(cls, text):
        lines = [line.strip() for line in str(text).splitlines() if line.strip()]
        if not lines:
            return cls
        headers = lines[0].split()
        if not headers or headers[0].lower() != "name":
            return cls
        for line in lines[1:]:
            values = line.split()
            if len(values) < 6:
                continue
            entry = _VirtualType2Record(*values[:6])
            cls._types[cls._norm(values[0])] = entry
        return cls


GlobalSetting.VirtualAtomTypes["vatom2"] = 3

__all__ = ["VirtualType2"]
