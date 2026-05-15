"""Legacy cmap_base compatibility helpers."""

from __future__ import annotations


class _CMapTypeRecord:
    def __init__(self, name, resolution, parameters):
        self.name = name
        self.resolution = int(resolution)
        self.parameters = list(parameters)


class CMapType:
    """Minimal executable compatibility for legacy CMAP type registry."""

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

    @staticmethod
    def same_force(atom_list):
        return [atom_list]

    @classmethod
    def New_From_Dict(cls, mapping):
        for name, payload in dict(mapping).items():
            if isinstance(payload, _CMapTypeRecord):
                entry = payload
            else:
                entry = _CMapTypeRecord(
                    name,
                    payload.get("resolution", 0),
                    payload.get("parameters", []),
                )
            cls._types[cls._norm(name)] = entry
        return cls


__all__ = ["CMapType"]
