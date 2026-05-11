"""RESP charge fitting using the same PySCF-based algorithm as Xponge."""

import math


default_radius = {
    "H": 1.2,
    "C": 1.5,
    "N": 1.5,
    "O": 1.4,
    "P": 1.8,
    "S": 1.75,
    "F": 1.35,
    "Cl": 1.7,
    "Br": 2.3,
}


def _require_numpy_pyscf():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError("NumPy is required for RESP charge calculation") from exc
    try:
        from pyscf import gto, scf
    except ImportError as exc:
        raise ImportError("PySCF is required for RESP charge calculation") from exc
    return np, gto, scf


def _fibonacci_grid(npoints, center, radius):
    if npoints <= 0:
        return []
    out = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(npoints):
        y = 1.0 - (2.0 * i + 1.0) / npoints
        r = math.sqrt(max(1.0 - y * y, 0.0))
        theta = golden_angle * i
        out.append(
            [
                center[0] + radius * math.cos(theta) * r,
                center[1] + radius * y,
                center[2] + radius * math.sin(theta) * r,
            ]
        )
    return out


def _get_mk_grid(assign, crd, area_density=1.0, layer=4, radius=None):
    np, _, _ = _require_numpy_pyscf()
    grids = []
    factor = area_density * 0.52918 * 0.52918 * 4 * np.pi
    real_radius = dict(default_radius)
    if radius:
        real_radius.update(radius)
    lists0 = np.array([1.4 + 0.2 * i for i in range(layer)])
    for i, atom in enumerate(assign.atoms):
        if atom not in real_radius:
            raise KeyError(f"Radius for element {atom} not found")
        r0 = real_radius[atom] / 0.52918
        for r in r0 * lists0:
            grids.extend(_fibonacci_grid(int(factor * r * r), crd[i], r))
    grids = np.array(grids).reshape(-1, 3)
    for i, atom in enumerate(assign.atoms):
        r0 = 1.39 * real_radius[atom] / 0.52918
        t = np.linalg.norm(grids - crd[i], axis=1)
        grids = grids[t >= r0, :]
    return grids


def _force_equivalence_q(q, extra_equivalence):
    np, _, _ = _require_numpy_pyscf()
    for eq_group in extra_equivalence:
        q_mean = np.mean(q[eq_group])
        q[eq_group] = q_mean
    return [float(x) for x in np.asarray(q).reshape(-1)]


def _get_pyscf_mol(assign, basis, charge, spin, opt):
    _, gto, scf = _require_numpy_pyscf()
    mols = ""
    for i, atom in enumerate(assign.atoms):
        x, y, z = assign.coordinates[i]
        mols += f"{atom} {x:f} {y:f} {z:f}\n"
    mol = gto.M(atom=mols, verbose=0, basis=basis, charge=charge, spin=spin)
    fun = scf.RHF(mol) if spin == 0 else scf.UHF(mol)
    if opt:
        from pyscf.geomopt.geometric_solver import optimize as geometric_opt

        mol = geometric_opt(fun)
        fun = scf.RHF(mol) if spin == 0 else scf.UHF(mol)
        for i, coord in enumerate(mol.atom_coords() * 0.52918):
            assign.set_coordinate(i, float(coord[0]), float(coord[1]), float(coord[2]))
    fun.run()
    return mol, fun


def _resp_scf_kernel(mol, assign, a, b, matrix_a, matrix_a0, matrix_b, q):
    np, _, _ = _require_numpy_pyscf()
    step = 0
    q_last_step = q
    while step == 0 or np.max(np.abs(q - q_last_step)) > 1e-4:
        step += 1
        q_last_step = q
        for i in range(mol.natm):
            if assign.atoms[i] != "H":
                matrix_a[i][i] = matrix_a0[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
        q = np.dot(np.linalg.inv(matrix_a), matrix_b)[:-1]
    return q


def _find_tofit_second(mol, assign):
    tofit_second = []
    fit_group = {i: -1 for i in range(mol.natm)}
    sublength = 0
    for i in range(mol.natm):
        if _atom_judge(assign, i, "C4"):
            fit_group[i] = len(tofit_second)
            tofit_second.append([i])
            hydrogens = [j for j in assign.bonds[i] if assign.atoms[j] == "H"]
            if hydrogens:
                for j in hydrogens:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(hydrogens)
                sublength += len(hydrogens) - 1
        if _atom_judge(assign, i, "C3"):
            hydrogens = [j for j in assign.bonds[i] if assign.atoms[j] == "H"]
            if len(hydrogens) == 2:
                fit_group[i] = len(tofit_second)
                tofit_second.append([i])
                for j in hydrogens:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(hydrogens)
                sublength += 1
    return tofit_second, fit_group, sublength


def _atom_judge(assign, atom, mask):
    element = "".join(ch for ch in mask if not ch.isdigit())
    digits = "".join(ch for ch in mask if ch.isdigit())
    if digits:
        return assign.atoms[atom] == element and len(assign.bonds[atom]) == int(digits)
    return assign.atoms[atom] == element


def _correct_extra_equivalence(tofit_second, fit_group, sublength, extra_equivalence, atom_numbers):
    if not extra_equivalence:
        return tofit_second, fit_group, sublength
    equi_group = []
    for eq in extra_equivalence:
        group = sorted({fit_group[eq_atom] for eq_atom in eq if fit_group[eq_atom] != -1})
        equi_group.append(group)
    all_groups_list = sorted({fit_group[atom] for atom in range(atom_numbers)})
    group_map = {i: i for i in all_groups_list}
    for eq in equi_group:
        for group in eq:
            group_map[group] = eq[0]
    temp_max = 0
    for group in all_groups_list:
        if group == -1:
            continue
        if group_map[group] == group:
            group_map[group] = temp_max
            temp_max += 1
        else:
            group_map[group] = group_map[group_map[group]]
    old = tofit_second
    tofit_second = [[] for _ in range(temp_max)]
    for i, group in enumerate(old):
        tofit_second[group_map[i]].extend(group)
        sublength -= len(group) - 1
    for group in tofit_second:
        sublength += len(group) - 1
    for atom in range(atom_numbers):
        fit_group[atom] = group_map[fit_group[atom]]
    return tofit_second, fit_group, sublength


def _get_a20_and_b20(total_length, tofit_second, fit_group, sublength, mol, matrix_a0, matrix_b, charge, q):
    np, _, _ = _require_numpy_pyscf()
    a20 = np.zeros((total_length, total_length))
    count = len(tofit_second)
    for i in range(mol.natm):
        if fit_group[i] == -1:
            fit_group[i] = count
            count += 1
        a20[mol.natm - sublength][fit_group[i]] += 1
        a20[fit_group[i]][mol.natm - sublength] += 1
    b20 = np.zeros(total_length)
    for i in range(mol.natm):
        b20[fit_group[i]] += matrix_b[i]
        for j in range(mol.natm):
            a20[fit_group[i]][fit_group[j]] += matrix_a0[i][j]
    b20[mol.natm - sublength] = charge
    count = 0
    for i in range(mol.natm):
        if fit_group[i] >= len(tofit_second):
            b20[mol.natm - sublength + count + 1] = q[i]
            a20[mol.natm - sublength + count + 1][len(tofit_second) + count] = 1
            a20[len(tofit_second) + count][mol.natm - sublength + count + 1] = 1
            count += 1
    return a20, b20


def resp_fit(assign, basis="6-31g*", opt=False, charge=None, spin=0, extra_equivalence=None,
             grid_density=6, grid_cell_layer=4, radius=None, a1=0.0005, a2=0.001,
             two_stage=True, only_esp=False):
    np, gto, _ = _require_numpy_pyscf()
    if extra_equivalence is None:
        extra_equivalence = []
    if charge is None:
        charge = int(round(sum(assign.charges)))

    mol, fun = _get_pyscf_mol(assign, basis, charge, spin, opt)
    grids = _get_mk_grid(assign, mol.atom_coords(), grid_density, grid_cell_layer, radius)
    vnuc = 0
    matrix_a0 = np.zeros((mol.natm, mol.natm))
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        z = mol.atom_charge(i)
        rp = r - grids
        for j in range(mol.natm):
            rpj = mol.atom_coord(j) - grids
            matrix_a0[i][j] = np.sum(1.0 / np.linalg.norm(rp, axis=1) / np.linalg.norm(rpj, axis=1))
        vnuc += z / np.einsum("xi,xi->x", rp, rp) ** 0.5

    matrix_a0 = np.hstack((matrix_a0, np.ones(mol.natm).reshape(-1, 1)))
    temp = np.ones(mol.natm + 1)
    temp[-1] = 0
    matrix_a0 = np.vstack((matrix_a0, temp.reshape(1, -1)))
    matrix_a = np.zeros_like(matrix_a0)
    matrix_a[:] = matrix_a0

    try:
        from pyscf import df

        fakemol = gto.fakemol_for_charges(grids)
        vele = np.einsum("ijp,ij->p", df.incore.aux_e2(mol, fakemol), fun.make_rdm1())
    except MemoryError:
        dm = fun.make_rdm1()
        vele = []
        for p in grids:
            mol.set_rinv_orig_(p)
            vele.append(np.einsum("ij,ij", mol.intor("int1e_rinv"), dm))
        vele = np.array(vele)
    mep = vnuc - vele

    matrix_b = np.zeros((mol.natm + 1))
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        rp = np.linalg.norm(r - grids, axis=1)
        matrix_b[i] = np.sum(mep / rp)
    matrix_b[-1] = charge
    matrix_b = matrix_b.reshape(-1, 1)
    q = np.dot(np.linalg.inv(matrix_a), matrix_b.reshape(-1, 1))[:-1]

    if only_esp:
        return _force_equivalence_q(q, extra_equivalence)

    q = _resp_scf_kernel(mol, assign, a1, 0.1, matrix_a, matrix_a0, matrix_b, q)
    if not two_stage:
        return _force_equivalence_q(q, extra_equivalence)

    tofit_second, fit_group, sublength = _find_tofit_second(mol, assign)
    tofit_second, fit_group, sublength = _correct_extra_equivalence(
        tofit_second, fit_group, sublength, extra_equivalence, mol.natm
    )
    if tofit_second:
        total_length = mol.natm - sublength + 1 + mol.natm - sublength - len(tofit_second)
        a20, b20 = _get_a20_and_b20(total_length, tofit_second, fit_group, sublength, mol, matrix_a0, matrix_b, charge, q)
        matrix_a = np.zeros_like(a20)
        matrix_a[:] = a20[:]
        matrix_b = b20.reshape(-1, 1)
        q_temp = np.dot(np.linalg.inv(matrix_a), matrix_b)[:-1]
        step = 0
        q_last_step = q_temp
        while step == 0 or np.max(np.abs(q_temp - q_last_step)) > 1e-4:
            step += 1
            q_last_step = q_temp
            for i in range(mol.natm - sublength):
                if assign.atoms[i] != "H":
                    matrix_a[i][i] = a20[i][i] + a2 / np.sqrt(q_last_step[i] * q_last_step[i] + 0.1 * 0.1)
            q_temp = np.dot(np.linalg.inv(matrix_a), matrix_b)[:-1]
        for i, group in enumerate(tofit_second):
            for j in group:
                q[j] = q_temp[i]
    return _force_equivalence_q(q, extra_equivalence)
