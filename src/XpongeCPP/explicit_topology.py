"""Validated explicit SPONGE topology subset, in angstrom/radian/kcal/mol.

Periodic torsions are ordered energy terms, including periodic impropers;
they are never re-inferred from a chemical graph. See docs/explicit-topology.md.
"""
import math
import operator
import hashlib


def _index(value, count):
    if isinstance(value, bool):
        raise ValueError("Boolean atom index")
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise ValueError("Expected an integer atom index") from exc
    if not 0 <= value < count:
        raise ValueError("Atom index out of range")
    return value


def _number(value):
    if isinstance(value, (bool, str)):
        raise ValueError("Expected a real parameter")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Nonfinite parameter")
    return value


def _validate(topology, count):
    widths = {"bonds": (2, 4), "angles": (3, 5),
              "dihedrals": (4, 7), "nb14s": (2, 4)}
    if set(topology) != set(widths) | {"exclusions"}:
        raise ValueError("Expected bonds, angles, dihedrals, nb14s and exclusions")
    result = {}
    for key, (natom, width) in widths.items():
        result[key] = []
        for row in topology[key]:
            if len(row) != width:
                raise ValueError("Invalid " + key + " row width")
            ids = [_index(v, count) for v in row[:natom]]
            if len(set(ids)) != natom:
                raise ValueError("Repeated atom in interaction")
            params = [_number(v) for v in row[natom:]]
            if key == "dihedrals":
                params[0] = _index(row[4], 2**31)
                if params[0] < 1:
                    raise ValueError("Nonpositive periodicity")
            elif params[0] < 0:
                raise ValueError("Negative force constant/scaling")
            if key == "bonds" and params[1] <= 0:
                raise ValueError("Nonpositive bond length")
            if key == "angles" and not 0 <= params[1] <= math.pi:
                raise ValueError("Angle outside [0, pi]")
            if key == "nb14s" and params[1] < 0:
                raise ValueError("Negative electrostatic scaling")
            result[key].append(ids + params)
    rows = topology["exclusions"]
    if len(rows) != count:
        raise ValueError("Exclusions must have one row per atom")
    exclusions = []
    for i, row in enumerate(rows):
        ids = [_index(v, count) for v in row]
        if i in ids or len(set(ids)) != len(ids):
            raise ValueError("Self or duplicate exclusion")
        exclusions.append(set(ids))
    for i, row in enumerate(exclusions):
        if any(i not in exclusions[j] for j in row):
            raise ValueError("Exclusions must be symmetric")
    result["exclusions"] = [sorted(row) for row in exclusions]
    return result


def _lj_identity(epsilon, rmin):
    epsilon, rmin = _number(epsilon), _number(rmin)
    if epsilon < 0 or rmin < 0:
        raise ValueError("Negative LJ parameter")
    name = "EXPLICIT_LJ_" + hashlib.sha256(
        (epsilon.hex() + ":" + rmin.hex()).encode()).hexdigest()
    return name, epsilon, rmin


def get_explicit_topology(molecule):
    """Return an owned copy of the five explicit topology tables."""
    return {k: [list(r) for r in rows]
            for k, rows in molecule._get_explicit_topology().items()}


def set_explicit_topology(molecule, topology):
    """Validate then atomically replace the five topology tables."""
    molecule._set_explicit_topology(_validate(topology, molecule.atom_count))


def set_atom_lj(molecule, atom_index, epsilon, rmin):
    """Set LJ epsilon (kcal/mol), Rmin/2 (angstrom), using a content-addressed type.

    This finalizes atom typing; call after topology materialization. New bonded
    parameters must be supplied explicitly if the graph is changed afterwards.
    """
    from . import register_amber_lj_parameter
    atom_index = _index(atom_index, molecule.atom_count)
    name, epsilon, rmin = _lj_identity(epsilon, rmin)
    if not molecule.has_topology_override:
        raise ValueError("Set explicit topology before assigning atom LJ")
    register_amber_lj_parameter(name, name, epsilon, rmin)
    molecule._set_lj_parameter_override(name, name, epsilon, rmin, "explicit-topology")
    molecule.atoms[atom_index].type = name
    return name
