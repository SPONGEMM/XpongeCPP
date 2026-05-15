"""Legacy angle_base compatibility helpers."""

from __future__ import annotations


class _AngleTypeRecord:
    def __init__(self, name, k, b):
        self.name = name
        self.k = float(k)
        self.b = float(b)


class AngleType:
    """Minimal executable compatibility for legacy AngleType usage."""

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
            cls._types[cls._norm(values[0])] = _AngleTypeRecord(values[0], values[1], values[2])
        return cls


__all__ = ["AngleType"]
