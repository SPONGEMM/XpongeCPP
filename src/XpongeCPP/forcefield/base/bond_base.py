"""Legacy bond_base compatibility helpers."""

from __future__ import annotations


class _BondTypeRecord:
    def __init__(self, name, k, b):
        self.name = name
        self.k = float(k)
        self.b = float(b)


class _BondForceEntity:
    def __init__(self, atoms, bond_type):
        self.atoms = list(atoms)
        self.type = bond_type
        self.k = bond_type.k
        self.b = bond_type.b


class BondType:
    """Minimal executable compatibility for legacy BondType usage."""

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
    def entity(cls, atoms, bond_type):
        return _BondForceEntity(atoms, bond_type)

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
            cls._types[cls._norm(values[0])] = _BondTypeRecord(values[0], values[1], values[2])
        return cls


def _gmx_parser(words, mol, stat):
    def _type_name(type_):
        return type_ if isinstance(type_, str) else type_.name

    atom1 = stat[int(words[0])]
    atom2 = stat[int(words[1])]
    if len(words) == 3:
        type1 = _type_name(atom1.type)
        type2 = _type_name(atom2.type)
        string = f"{type1}-{type2}"
        reversed_string = f"{type2}-{type1}"
        if BondType._norm(string) in BondType.get_all_types():
            mol.add_bonded_force(BondType.entity([atom1, atom2], BondType.getType(string)))
        elif BondType._norm(reversed_string) in BondType.get_all_types():
            mol.add_bonded_force(BondType.entity([atom1, atom2], BondType.getType(reversed_string)))
        else:
            raise KeyError(
                "Bonded Force Type "
                f"{string} (or {reversed_string}) not found. "
                "Did you import the proper force field / load bondtypes, "
                "or provide explicit bond parameters (5 columns) in the [ bonds ] section?"
            )
    elif len(words) == 5:
        new_force = BondType.entity([atom1, atom2], _BondTypeRecord("UNKNOWNS", float(words[3]), float(words[4])))
        mol.add_bonded_force(new_force)
    else:
        raise ValueError(f"Only 3 or 5 words should be in the line '{' '.join(words)}'")


__all__ = ["BondType", "_gmx_parser"]
