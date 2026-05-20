"""Legacy helper shim for XpongeCPP."""

from __future__ import annotations

from abc import ABC

from .._compat.imports import (
    Generate_New_Bonded_Force_Type,
    Generate_New_Pairwise_Force_Type,
    Xdict,
    Xopen,
    Xpri,
    Xprint,
    debug,
    generate_new_bonded_force_type,
    generate_new_pairwise_force_type,
    remove_real_global_variable,
    source,
    set_global_alternative_names,
    set_real_global_variable,
    xopen,
    xprint,
)
from .._core import Atom, Molecule, Residue, ResidueType
from ..gromacs import GlobalSetting
from .cv import CVSystem
from .file import file_filter, import_python_script, pdb_filter
from .math import (
    Guess_Element_From_Mass,
    Get_Basis_Vectors_From_Length_And_Angle,
    Get_Fibonacci_Grid,
    Get_Length_Angle_From_Basis_Vectors,
    Get_Rotate_Matrix,
    Kabsch,
    get_basis_vectors_from_length_and_angle,
    get_fibonacci_grid,
    get_length_angle_from_basis_vectors,
    get_rotate_matrix,
    guess_element_from_mass,
    kabsch,
)


class Type(ABC):
    """Minimal legacy-compatible base type class."""

    @classmethod
    def get_class_name(cls):
        return cls.__name__.replace("Type", "") or cls.__name__

    def __repr__(self):
        name = getattr(self, "name", self.__class__.__name__)
        return f"Type of {self.get_class_name()}: {name}"


class AbstractMolecule(ABC):
    """Minimal legacy-compatible abstract molecule marker."""


class Entity(ABC):
    """Minimal legacy-compatible entity marker."""


class ResidueLink:
    """Minimal legacy-compatible residue-link object."""

    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.tohash = ResidueLink.get_hash(atom1, atom2)

    def __repr__(self):
        return "Entity of ResidueLink: " + repr(self.atom1) + "-" + repr(self.atom2)

    def __hash__(self):
        return hash(self.tohash)

    def __eq__(self, other):
        return isinstance(other, ResidueLink) and self.tohash == other.tohash

    @staticmethod
    def get_hash(atom1, atom2):
        ids = sorted((id(atom1), id(atom2)))
        return tuple(ids)

    def deepcopy(self, forcopy):
        del forcopy
        return ResidueLink(self.atom1, self.atom2)


class AtomType(Type):
    """Minimal legacy-compatible AtomType registry placeholder."""

    _types = Xdict(not_found_message="Atom type {} not found. Did you import the proper force field?")
    _parameters = {"name": str}
    _property_units = {}

    def __init__(self, name, **kwargs):
        self.name = str(name)
        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def get_type(cls, name):
        return cls._types[str(name)]

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
    def Add_Property(cls, extra_properties):
        cls._parameters.update(dict(extra_properties))

    @classmethod
    def Set_Property_Unit(cls, property_name, category, unit):
        cls._property_units[str(property_name)] = (str(category), str(unit))

    @classmethod
    def Get_Property_Unit(cls, property_name):
        return cls._property_units.get(str(property_name))

    @classmethod
    def New_From_String(cls, text):
        lines = [line.strip() for line in str(text).splitlines() if line.strip()]
        if not lines:
            return cls
        entries = []
        if len(lines) > 1 and lines[0].split()[0].lower() == "name":
            headers = lines[0].split()
            cls._parameters.update({header: str for header in headers})
            for line in lines[1:]:
                values = line.split()
                if not values:
                    continue
                mapping = {}
                for index, header in enumerate(headers[1:], start=1):
                    if index < len(values):
                        raw = values[index]
                        try:
                            mapping[header] = float(raw)
                        except ValueError:
                            mapping[header] = raw
                entries.append((values[0], mapping))
        else:
            for line in lines:
                for name in line.split():
                    entries.append((name, {}))
        for name, mapping in entries:
            cls._types[name] = cls(name, **mapping)
        return cls


Entity.register(Atom)
Entity.register(Residue)
Entity.register(ResidueLink)
Entity.register(Molecule)
AbstractMolecule.register(Residue)
AbstractMolecule.register(ResidueType)
AbstractMolecule.register(Molecule)

__all__ = [
    "AbstractMolecule",
    "Atom",
    "AtomType",
    "CVSystem",
    "Entity",
    "Generate_New_Bonded_Force_Type",
    "Generate_New_Pairwise_Force_Type",
    "GlobalSetting",
    "Molecule",
    "Residue",
    "ResidueLink",
    "ResidueType",
    "Type",
    "Xdict",
    "Xopen",
    "Xpri",
    "Xprint",
    "debug",
    "file_filter",
    "get_basis_vectors_from_length_and_angle",
    "get_fibonacci_grid",
    "get_length_angle_from_basis_vectors",
    "get_rotate_matrix",
    "guess_element_from_mass",
    "import_python_script",
    "kabsch",
    "pdb_filter",
    "Guess_Element_From_Mass",
    "Get_Basis_Vectors_From_Length_And_Angle",
    "Get_Fibonacci_Grid",
    "Get_Length_Angle_From_Basis_Vectors",
    "Get_Rotate_Matrix",
    "Kabsch",
    "generate_new_bonded_force_type",
    "generate_new_pairwise_force_type",
    "remove_real_global_variable",
    "source",
    "set_global_alternative_names",
    "set_real_global_variable",
    "xopen",
    "xprint",
]
