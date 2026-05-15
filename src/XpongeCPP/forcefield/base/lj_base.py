"""Legacy lj_base compatibility helpers."""

from __future__ import annotations

import math


class _LJEntry:
    def __init__(self, name, epsilon, rmin):
        self.name = name
        self.epsilon = float(epsilon)
        self.rmin = float(rmin)


class LJType:
    """Minimal executable compatibility for legacy LJType lookups."""

    combining_method_A = None
    combining_method_B = None
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
    def _store(cls, name, epsilon, rmin):
        entry = _LJEntry(name, epsilon, rmin)
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
        keyset = {header.lower() for header in headers}
        for line in lines[1:]:
            values = line.split()
            if not values:
                continue
            row = {headers[i].lower(): values[i] for i in range(min(len(headers), len(values)))}
            name = row["name"]
            if {"a", "b"} <= keyset:
                a = float(row["a"])
                b = float(row["b"])
                if a == 0.0 or b == 0.0:
                    sigma = 0.0
                    epsilon = 0.0
                else:
                    sigma = (a / b) ** (1.0 / 6.0)
                    epsilon = 0.25 * b * sigma ** (-6.0)
                rmin = sigma * (4 ** (1 / 12) / 2)
            elif {"epsilon", "rmin"} <= keyset:
                epsilon = float(row["epsilon"])
                rmin = float(row["rmin"])
            elif {"epsilon[ev]", "sigma[nm]"} <= keyset:
                epsilon = float(row["epsilon[ev]"]) * 23.06054783061903
                sigma = float(row["sigma[nm]"]) * 10.0
                rmin = sigma * (4 ** (1 / 12) / 2)
            else:
                continue
            cls._store(name, epsilon, rmin)
        return cls


def Lorentz_Berthelot_For_A(epsilon1, rmin1, epsilon2, rmin2):
    return math.sqrt(epsilon1 * epsilon2) * ((rmin1 + rmin2) ** 12)


def Lorentz_Berthelot_For_B(epsilon1, rmin1, epsilon2, rmin2):
    return math.sqrt(epsilon1 * epsilon2) * 2 * ((rmin1 + rmin2) ** 6)


lorentz_berthelot_for_a = Lorentz_Berthelot_For_A
lorentz_berthelot_for_b = Lorentz_Berthelot_For_B

__all__ = [
    "LJType",
    "Lorentz_Berthelot_For_A",
    "Lorentz_Berthelot_For_B",
    "lorentz_berthelot_for_a",
    "lorentz_berthelot_for_b",
]
