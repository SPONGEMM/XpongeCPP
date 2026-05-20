"""PySCF RESP backend helpers.

This module preserves the current QM behavior while the RESP migration is in
progress. It exists so the orchestration layer can switch to Psi4 later without
changing the Python compatibility API.
"""

from __future__ import annotations

import time


def require_numpy_pyscf():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError("NumPy is required for RESP charge calculation") from exc
    try:
        from pyscf import gto, scf
    except ImportError as exc:
        raise ImportError("PySCF is required for RESP charge calculation") from exc
    return np, gto, scf


def _build_pyscf_atom_block(assign):
    lines = []
    for i, atom in enumerate(assign.atoms):
        x, y, z = assign.coordinates[i]
        lines.append(f"{atom} {x:f} {y:f} {z:f}")
    return "\n".join(lines)


def _build_pyscf_molecule(assign, basis, charge, spin, gto):
    return gto.M(atom=_build_pyscf_atom_block(assign), verbose=0, basis=basis, charge=charge, spin=spin)


def _build_pyscf_wavefunction(mol, spin, scf):
    return scf.RHF(mol) if spin == 0 else scf.UHF(mol)


def build_backend_payload(assign, basis, charge, spin, opt, return_timings=False):
    np, gto, scf = require_numpy_pyscf()
    timings = {}
    total_start = time.perf_counter()
    start = time.perf_counter()
    mol = _build_pyscf_molecule(assign, basis, charge, spin, gto)
    timings["build"] = time.perf_counter() - start
    start = time.perf_counter()
    wavefunction = _build_pyscf_wavefunction(mol, spin, scf)
    if opt:
        from pyscf.geomopt.geometric_solver import optimize as geometric_opt

        mol = geometric_opt(wavefunction)
        wavefunction = _build_pyscf_wavefunction(mol, spin, scf)
        for i, coord in enumerate(mol.atom_coords() * 0.52918):
            assign.set_coordinate(i, float(coord[0]), float(coord[1]), float(coord[2]))
    wavefunction.run()
    timings["scf"] = time.perf_counter() - start
    atom_coordinates_bohr = np.array([mol.atom_coord(i) for i in range(mol.natm)])
    nuclear_charges = np.array([mol.atom_charge(i) for i in range(mol.natm)], dtype=float)
    payload = {
        "atom_symbols": list(assign.atoms),
        "atom_coordinates_bohr": atom_coordinates_bohr,
        "nuclear_charges": nuclear_charges,
        "mol": mol,
        "wavefunction": wavefunction,
        "charge": charge,
        "spin": spin,
    }
    if return_timings:
        timings["total"] = time.perf_counter() - total_start
        payload["timings"] = timings
    return payload


def compute_esp_on_grid(payload, grid_points_bohr):
    np, gto, _ = require_numpy_pyscf()
    grid_points_bohr = np.asarray(grid_points_bohr, dtype=float)
    mol = payload["mol"]
    wavefunction = payload["wavefunction"]
    try:
        from pyscf import df

        fakemol = gto.fakemol_for_charges(grid_points_bohr)
        return np.einsum("ijp,ij->p", df.incore.aux_e2(mol, fakemol), wavefunction.make_rdm1())
    except MemoryError:
        dm = wavefunction.make_rdm1()
        vele = []
        for point in grid_points_bohr:
            mol.set_rinv_orig_(point)
            vele.append(np.einsum("ij,ij", mol.intor("int1e_rinv"), dm))
        return np.array(vele)
