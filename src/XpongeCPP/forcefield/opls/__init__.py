"""OPLS force-field compatibility helpers."""

import math
from pathlib import Path

from ... import (
    clear_amber_dihedral_parameters,
    clear_amber_improper_parameters,
    register_amber_angle_parameter,
    register_amber_bond_parameter,
    register_amber_improper_dihedral_parameter,
    register_amber_lj_parameter,
    register_amber_nb14_scale,
    register_amber_proper_dihedral_parameter,
    set_lj_combining_rule,
)
from ...gromacs import load_ffitp
from .. import repository_reference_forcefield_path


OPLS_BOND_TYPE_MAP = {}
_KJ_PER_KCAL = 4.184
_NM_PER_ANGSTROM = 0.1


def data_path(*parts):
    return repository_reference_forcefield_path("opls", *parts)


def _iter_parameter_rows(text):
    for index, line in enumerate(text.splitlines()):
        if index == 0:
            continue
        line = line.strip()
        if line:
            yield line.split()


def _split_atom_name(name, expected):
    parts = [part.strip() for part in name.split("-")]
    if len(parts) != expected:
        raise ValueError(f"expected {expected} atom labels in {name!r}")
    return parts


def _sigma_to_rmin(sigma):
    return sigma * (4 ** (1 / 12) / 2)


def _kj_to_kcal(value):
    return value / _KJ_PER_KCAL


def _nm_to_angstrom(value):
    return value / _NM_PER_ANGSTROM


def _bond_k_to_internal(value):
    return value / (_KJ_PER_KCAL * (_NM_PER_ANGSTROM ** -2))


def _degrees_to_radians(value):
    return math.radians(value)


def load_parameter_from_ffitp(filename, folder, reset=True):
    filename = Path(folder) / filename
    output = load_ffitp(str(filename))
    set_lj_combining_rule("good_hope")
    if reset:
        clear_amber_dihedral_parameters()
        clear_amber_improper_parameters()

    OPLS_BOND_TYPE_MAP.clear()
    OPLS_BOND_TYPE_MAP.update(output.get("bond_type_names", {}))

    for words in _iter_parameter_rows(output["LJ"]):
        atom_type1, atom_type2 = _split_atom_name(words[0], 2)
        if atom_type1 != atom_type2:
            continue
        register_amber_lj_parameter(
            atom_type1,
            atom_type1,
            _kj_to_kcal(float(words[2])),
            _nm_to_angstrom(_sigma_to_rmin(float(words[1]))),
        )

    for words in _iter_parameter_rows(output["bonds"]):
        register_amber_bond_parameter(
            *_split_atom_name(words[0], 2),
            _bond_k_to_internal(float(words[2])),
            _nm_to_angstrom(float(words[1])),
        )

    for words in _iter_parameter_rows(output["angles"]):
        register_amber_angle_parameter(
            _split_atom_name(words[0], 3),
            _kj_to_kcal(float(words[2])),
            _degrees_to_radians(float(words[1])),
        )

    for words in _iter_parameter_rows(output["dihedrals"]):
        register_amber_proper_dihedral_parameter(
            _split_atom_name(words[0], 4),
            int(words[3]),
            _kj_to_kcal(float(words[2])),
            _degrees_to_radians(float(words[1])),
            bool(int(words[4])),
        )

    for words in _iter_parameter_rows(output["periodic_impropers"]):
        register_amber_improper_dihedral_parameter(
            _split_atom_name(words[0], 4),
            int(words[3]),
            _kj_to_kcal(float(words[2])),
            _degrees_to_radians(float(words[1])),
        )

    for words in _iter_parameter_rows(output["nb14"]):
        register_amber_nb14_scale(*_split_atom_name(words[0], 2), float(words[1]), float(words[2]))

    return output


Load_Parameter_From_FFITP = load_parameter_from_ffitp
