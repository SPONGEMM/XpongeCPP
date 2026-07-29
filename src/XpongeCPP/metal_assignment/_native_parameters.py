"""Extract assigned parameters from the native SPONGE raw-export surface."""

from __future__ import annotations

import math
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Sequence


def _lines(path: Path) -> list[list[str]]:
    return [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _canonical_atoms(kind: str, atom_ids: Sequence[str]) -> tuple[str, ...]:
    values = tuple(atom_ids)
    if kind in {"bond", "angle", "proper_dihedral"}:
        return min(values, values[::-1])
    center = values[2]
    outers = sorted((values[0], values[1], values[3]))
    return (outers[0], center, outers[1], outers[2])


def _lj_parameters(path: Path, atom_ids: Sequence[str], source: str):
    tokens = path.read_text(encoding="utf-8").split()
    atom_count = int(tokens[0])
    type_count = int(tokens[1])
    if atom_count != len(atom_ids):
        raise ValueError("native LJ export atom coverage mismatch")
    triangle_count = type_count * (type_count + 1) // 2
    cursor = 2
    coefficients_a = [float(value) for value in tokens[cursor:cursor + triangle_count]]
    cursor += triangle_count
    coefficients_b = [float(value) for value in tokens[cursor:cursor + triangle_count]]
    cursor += triangle_count
    assignments = [int(value) for value in tokens[cursor:cursor + atom_count]]
    if len(assignments) != atom_count:
        raise ValueError("native LJ export type coverage mismatch")
    by_type = {}
    for type_index in range(type_count):
        diagonal = type_index * (type_index + 1) // 2 + type_index
        coefficient_a = coefficients_a[diagonal]
        coefficient_b = coefficients_b[diagonal]
        if coefficient_a == 0.0 or coefficient_b == 0.0:
            epsilon = 0.0
            rmin = 0.0
        else:
            sigma = (coefficient_a / coefficient_b) ** (1.0 / 6.0)
            epsilon = 0.25 * coefficient_b * sigma ** -6.0
            rmin = sigma * (4.0 ** (1.0 / 12.0) / 2.0)
        if not all(math.isfinite(value) for value in (epsilon, rmin)):
            raise ValueError("native LJ export contains non-finite parameters")
        by_type[type_index] = {
            "epsilon": epsilon,
            "rmin": rmin,
            "energy_unit": "kcal/mol",
            "length_unit": "angstrom",
            "source": source,
        }
    return {
        atom_id: dict(by_type[assignments[index]])
        for index, atom_id in enumerate(atom_ids)
    }


def extract_native_parameters(restype: Any, atom_ids: Sequence[str], source: str):
    """Return per-atom LJ values and neutral bonded-parameter records."""
    from XpongeCPP import (
        Molecule,
        get_template_molecule,
        has_template,
        molecule_from_residuetype,
        save_sponge_input_raw,
    )

    with TemporaryDirectory(prefix="xponge-native-parameters-") as tempdir:
        prefix = Path(tempdir) / "assigned"
        if (
            isinstance(restype, Molecule)
            or (
                hasattr(restype, "residue_count")
                and hasattr(restype, "residues")
                and hasattr(restype, "atoms")
            )
        ):
            molecule = restype
        elif hasattr(restype, "name") and has_template(restype.name):
            molecule = get_template_molecule(restype.name)
        else:
            molecule = molecule_from_residuetype(restype)
        save_sponge_input_raw(molecule, str(prefix))
        atom_types = [
            fields[0]
            for fields in _lines(Path(f"{prefix}_atom_type_name.txt"))[1:]
        ]
        if len(atom_types) != len(atom_ids):
            raise ValueError("native atom-type export coverage mismatch")
        lj_parameters = _lj_parameters(
            Path(f"{prefix}_LJ.txt"), atom_ids, source
        )

        records: list[dict[str, Any]] = []
        connectivity: set[tuple[int, int]] = set()
        bond_lines = _lines(Path(f"{prefix}_bond.txt"))
        for fields in bond_lines[1:]:
            atom1, atom2 = int(fields[0]), int(fields[1])
            connectivity.add(tuple(sorted((atom1, atom2))))
            record_atoms = _canonical_atoms(
                "bond", (atom_ids[atom1], atom_ids[atom2])
            )
            records.append({
                "kind": "bond",
                "atom_ids": list(record_atoms),
                "parameters": {
                    "force_constant": float(fields[2]),
                    "equilibrium_distance": float(fields[3]),
                    "force_constant_unit": "kcal/mol/angstrom^2",
                    "distance_unit": "angstrom",
                    "type_name": f"{atom_types[atom1]}-{atom_types[atom2]}",
                },
                "source": source,
            })

        angle_lines = _lines(Path(f"{prefix}_angle.txt"))
        for fields in angle_lines[1:]:
            indices = tuple(int(value) for value in fields[:3])
            record_atoms = _canonical_atoms(
                "angle", tuple(atom_ids[index] for index in indices)
            )
            records.append({
                "kind": "angle",
                "atom_ids": list(record_atoms),
                "parameters": {
                    "force_constant": float(fields[3]),
                    "equilibrium_angle": float(fields[4]),
                    "force_constant_unit": "kcal/mol/radian^2",
                    "angle_unit": "radian",
                    "type_name": "-".join(atom_types[index] for index in indices),
                },
                "source": source,
            })

        grouped_dihedrals: dict[
            tuple[str, tuple[int, int, int, int]], list[tuple[int, float, float]]
        ] = {}
        dihedral_lines = _lines(Path(f"{prefix}_dihedral.txt"))
        for fields in dihedral_lines[1:]:
            indices = tuple(int(value) for value in fields[:4])
            is_proper = all(
                tuple(sorted(pair)) in connectivity
                for pair in zip(indices, indices[1:])
            )
            kind = "proper_dihedral" if is_proper else "improper_dihedral"
            grouped_dihedrals.setdefault((kind, indices), []).append(
                (int(fields[4]), float(fields[5]), float(fields[6]))
            )
        for (kind, indices), terms in grouped_dihedrals.items():
            record_atoms = _canonical_atoms(
                kind, tuple(atom_ids[index] for index in indices)
            )
            type_name = "-".join(atom_types[index] for index in indices)
            if kind == "proper_dihedral":
                parameters = {
                    "force_constants": [term[1] for term in terms],
                    "phases": [term[2] for term in terms],
                    "periodicities": [term[0] for term in terms],
                    "energy_unit": "kcal/mol",
                    "angle_unit": "radian",
                    "type_name": type_name,
                }
            else:
                periodicity, force_constant, phase = terms[0]
                parameters = {
                    "force_constant": force_constant,
                    "phase": phase,
                    "periodicity": periodicity,
                    "energy_unit": "kcal/mol",
                    "angle_unit": "radian",
                    "type_name": type_name,
                }
            records.append({
                "kind": kind,
                "atom_ids": list(record_atoms),
                "parameters": parameters,
                "source": source,
            })
    records.sort(
        key=lambda item: (
            item["kind"],
            item["atom_ids"],
            repr(sorted(item["parameters"].items())),
        )
    )
    return lj_parameters, records


__all__ = ["extract_native_parameters"]
