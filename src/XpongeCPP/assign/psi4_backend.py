"""Psi4 RESP backend helpers."""

from __future__ import annotations

import time

ANGSTROM_PER_BOHR = 0.52918


def require_numpy_psi4():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError("NumPy is required for RESP charge calculation") from exc
    try:
        import psi4
    except ImportError as exc:
        raise ImportError("Psi4 is required for RESP charge calculation") from exc
    return np, psi4


def _build_psi4_geometry_block(assign, charge, spin):
    multiplicity = int(spin) + 1
    geometry_lines = [f"{int(charge)} {multiplicity}"]
    for atom_index, atom in enumerate(assign.atoms):
        x, y, z = assign.coordinates[atom_index]
        geometry_lines.append(f"{atom} {x:.16f} {y:.16f} {z:.16f}")
    return "\n".join(geometry_lines)


def build_backend_payload(assign, basis, charge, spin, opt, return_timings=False):
    np, psi4 = require_numpy_psi4()
    timings = {}
    total_start = time.perf_counter()
    psi4.core.be_quiet()
    start = time.perf_counter()
    mol = psi4.geometry(_build_psi4_geometry_block(assign, charge, spin))
    psi4.core.clean_options()
    psi4.set_options(
        {
            "basis": basis,
            "reference": "rhf" if int(spin) == 0 else "uhf",
        }
    )
    timings["build"] = time.perf_counter() - start
    start = time.perf_counter()
    if opt:
        _, wavefunction = psi4.optimize("scf", molecule=mol, return_wfn=True)
    else:
        _, wavefunction = psi4.energy("scf", molecule=mol, return_wfn=True)
    timings["scf"] = time.perf_counter() - start
    out_mol = wavefunction.molecule()
    atom_coordinates_bohr = np.array(out_mol.geometry().np, dtype=float)
    nuclear_charges = np.array([out_mol.Z(i) for i in range(out_mol.natom())], dtype=float)
    if opt:
        for i, coord in enumerate(atom_coordinates_bohr * ANGSTROM_PER_BOHR):
            assign.set_coordinate(i, float(coord[0]), float(coord[1]), float(coord[2]))
    payload = {
        "atom_symbols": list(assign.atoms),
        "atom_coordinates_bohr": atom_coordinates_bohr,
        "nuclear_charges": nuclear_charges,
        "wavefunction": wavefunction,
        "charge": charge,
        "spin": spin,
    }
    if return_timings:
        timings["total"] = time.perf_counter() - total_start
        payload["timings"] = timings
    return payload


def compute_esp_on_grid(payload, grid_points_bohr):
    np, psi4 = require_numpy_psi4()
    wavefunction = payload["wavefunction"]
    grid_points_bohr = np.asarray(grid_points_bohr, dtype=float)
    grid_matrix = psi4.core.Matrix.from_array(grid_points_bohr * ANGSTROM_PER_BOHR)
    prop = psi4.core.ESPPropCalc(wavefunction)
    total_esp = np.array(prop.compute_esp_over_grid_in_memory(grid_matrix), dtype=float).reshape(-1)
    vnuc = np.zeros(len(grid_points_bohr), dtype=float)
    for coord, charge in zip(payload["atom_coordinates_bohr"], payload["nuclear_charges"]):
        rp = grid_points_bohr - coord
        vnuc += charge / np.linalg.norm(rp, axis=1)
    return vnuc - total_esp
