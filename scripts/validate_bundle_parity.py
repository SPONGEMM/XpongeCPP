#!/usr/bin/env python3
"""Generate and compare equivalent Xponge and XpongeCPP input bundles."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np


SYSTEM_BUILDER = r"""
import sys
from pathlib import Path

backend = sys.argv[1]
output_dir = Path(sys.argv[2])
xponge_root = Path(sys.argv[3])
if backend == "xponge":
    # XpongeCPP's editable install registers a meta-path finder for the legacy
    # package name. Remove only that finder so this process loads the reference
    # source tree while retaining the same binary Python dependencies.
    sys.meta_path = [
        finder for finder in sys.meta_path
        if "ScikitBuildRedirectingFinder" not in type(finder).__name__
    ]
    sys.path.insert(0, str(xponge_root))
    import Xponge as X
    import Xponge.forcefield.amber.ff14sb  # noqa: F401
else:
    import XpongeCPP as X
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

molecule = X.get_peptide_from_sequence("AA")
for index, atom in enumerate(molecule.atoms):
    atom.x = index * 0.1
    atom.y = index * 0.2
    atom.z = index * 0.3
molecule.box_length = [40.0, 41.0, 42.0]
if backend == "xponge":
    @X.Molecule.Set_Save_SPONGE_Input("bond_soft")
    def write_soft_bond(_):
        return "1\n0 1 12.5 1.25 1"

    @X.Molecule.Set_Save_SPONGE_Input("listed_forces")
    def write_listed_forces(_):
        return (
            "[[[ Ryckaert_Bellemans ]]]\n"
            "[[ parameters ]]\n"
            "int atom_a, int atom_b, int atom_c, int atom_d, "
            "float c0, float c1, float c2, float c3, float c4, float c5\n"
            "[[ potential ]]\n"
            "SADfloat<15> cphi = cosf(phi_abcd - CONSTANT_Pi);\n"
            "SADfloat<15> cphi2 = cphi * cphi;\n"
            "SADfloat<15> cphi3 = cphi2 * cphi;\n"
            "SADfloat<15> cphi4 = cphi3 * cphi;\n"
            "SADfloat<15> cphi5 = cphi4 * cphi;\n"
            "E = c0 + c1 * cphi + c2 * cphi2 + c3 * cphi3 + c4 * cphi4 + c5 * cphi5;\n"
            "[[ end ]]\n"
        )

    @X.Molecule.Set_Save_SPONGE_Input("Ryckaert_Bellemans")
    def write_ryckaert_bellemans(_):
        return "1\n1 2 3 4 0.1 0.2 0.3 0.4 0.5 0.6"
else:
    molecule.add_bond_soft(1, 0, 12.5, 1.25, 1)
    molecule.add_ryckaert_bellemans(4, 3, 2, 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
X.save_sponge_input_bundle(molecule, "system", output_dir)
"""


def _run_builder(backend: str, output_dir: Path, xponge_root: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-c", SYSTEM_BUILDER, backend, str(output_dir), str(xponge_root)],
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"{backend} bundle generation failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def _assert_close(lhs, rhs, label: str, *, rtol: float = 1e-5, atol: float = 1e-5) -> float:
    lhs = np.asarray(lhs)
    rhs = np.asarray(rhs)
    if lhs.shape != rhs.shape or not np.allclose(
        lhs, rhs, rtol=rtol, atol=atol, equal_nan=True
    ):
        raise AssertionError(f"{label} differs: shapes {lhs.shape} and {rhs.shape}")
    return float(np.max(np.abs(lhs - rhs))) if lhs.size else 0.0


def _assert_exact(lhs, rhs, label: str) -> None:
    lhs = np.asarray(lhs)
    rhs = np.asarray(rhs)
    if lhs.shape != rhs.shape or not np.array_equal(lhs, rhs):
        raise AssertionError(f"{label} differs: shapes {lhs.shape} and {rhs.shape}")


def _dataset_paths(handle: h5py.File) -> set[str]:
    paths: set[str] = set()
    handle.visititems(lambda name, obj: paths.add("/" + name) if isinstance(obj, h5py.Dataset) else None)
    return paths


def _content_hash(handle: h5py.File, bundle_file: str, path_prefixes=()) -> str:
    excluded = {
        "/schema/name",
        "/schema/version",
        "/parameters/sponge/schema/name",
        "/parameters/sponge/schema/version",
        "/topology/atom_count",
        "/topology/atom_order_hash",
        "/topology/topology_hash",
        "/topology/forcefield_hash",
        "/protocol/topology_compatibility/topology_hash",
        "/identity/content_hash",
        "/identity/uuid",
        "/protocol/cv_count",
        "/protocol/restraint_count",
    }
    datasets = {}
    handle.visititems(
        lambda name, item: datasets.__setitem__("/" + name, item)
        if isinstance(item, h5py.Dataset) and "/" + name not in excluded
        else None
    )
    digest = hashlib.sha256(bundle_file.encode())
    for path, dataset in sorted(datasets.items()):
        if path_prefixes and not path.startswith(path_prefixes):
            continue
        array = np.asarray(dataset[()])
        digest.update(b"\0" + path.encode() + b"\0")
        digest.update(str(dataset.dtype).encode("ascii") + b"\0")
        digest.update(repr(array.shape).encode("ascii") + b"\0")
        if array.dtype.kind in {"O", "U", "S"}:
            strings = np.asarray(dataset.asstr()[()]).reshape(-1)
            digest.update("\0".join(map(str, strings)).encode())
        else:
            digest.update(np.ascontiguousarray(array).tobytes())
    return "sha256:" + digest.hexdigest()


def _records(handle: h5py.File, group: str, parameters: tuple[str, ...]):
    atoms = handle[f"{group}/atoms"][:]
    values = [handle[f"{group}/{name}"][:] for name in parameters]
    records = []
    for index, row in enumerate(atoms):
        parameter_values = []
        for column in values:
            parameter_values.extend(float(value) for value in np.asarray(column[index]).reshape(-1))
        records.append((tuple(int(value) for value in row), tuple(parameter_values)))
    return sorted(
        records,
        key=lambda record: record[0] + tuple(round(value, 4) for value in record[1]),
    )


def _compare_records(lhs: h5py.File, rhs: h5py.File, group: str, parameters: tuple[str, ...]) -> float:
    lhs_records = _records(lhs, group, parameters)
    rhs_records = _records(rhs, group, parameters)
    if len(lhs_records) != len(rhs_records):
        raise AssertionError(f"{group} count differs")
    _assert_exact([record[0] for record in lhs_records], [record[0] for record in rhs_records], group)
    return _assert_close(
        [record[1] for record in lhs_records],
        [record[1] for record in rhs_records],
        f"{group} parameters",
    )


def _dihedral_records(handle: h5py.File):
    bonds = {tuple(sorted(map(int, row))) for row in handle["/forcefield/bond/atoms"][:]}
    atoms = handle["/forcefield/dihedral/atoms"][:]
    periodicity = handle["/forcefield/dihedral/periodicity"][:]
    force_constant = handle["/forcefield/dihedral/k"][:]
    phase = handle["/forcefield/dihedral/phi0"][:]
    records = []
    for index, row_values in enumerate(atoms):
        row = tuple(map(int, row_values))
        centers = [
            atom
            for atom in row
            if all(other == atom or tuple(sorted((atom, other))) in bonds for other in row)
        ]
        if len(centers) == 1:
            center = centers[0]
            atom_key = (1, center, *sorted(atom for atom in row if atom != center))
        else:
            atom_key = (0, *min(row, row[::-1]))
        parameters = (float(periodicity[index]), float(force_constant[index]), float(phase[index]))
        records.append((atom_key, parameters))
    return sorted(records, key=lambda record: record[0] + tuple(round(value, 4) for value in record[1]))


def _compare_dihedrals(lhs: h5py.File, rhs: h5py.File) -> float:
    lhs_records = _dihedral_records(lhs)
    rhs_records = _dihedral_records(rhs)
    _assert_exact([record[0] for record in lhs_records], [record[0] for record in rhs_records], "dihedrals")
    return _assert_close(
        [record[1] for record in lhs_records],
        [record[1] for record in rhs_records],
        "dihedral parameters",
    )


def _expanded_lj(handle: h5py.File) -> np.ndarray:
    atom_types = handle["/forcefield/lj/type"][:]
    pair_a = handle["/forcefield/lj/pair_A_12"][:]
    pair_b = handle["/forcefield/lj/pair_B_6"][:]

    def pair_value(values, lhs, rhs):
        high, low = max(int(lhs), int(rhs)), min(int(lhs), int(rhs))
        return values[high * (high + 1) // 2 + low]

    return np.asarray(
        [
            [
                (pair_value(pair_a, lhs, rhs), pair_value(pair_b, lhs, rhs))
                for rhs in atom_types
            ]
            for lhs in atom_types
        ]
    )


def _exclusion_pairs(handle: h5py.File) -> set[tuple[int, int]]:
    offsets = handle["/topology/exclusions/offset"][:]
    values = handle["/topology/exclusions/list"][:]
    return {
        tuple(sorted((atom, int(excluded))))
        for atom in range(len(offsets) - 1)
        for excluded in values[offsets[atom] : offsets[atom + 1]]
    }


def compare_bundles(reference_dir: Path, candidate_dir: Path) -> dict[str, object]:
    topology_name = "system_topology.spgt.h5"
    protocol_name = "system_protocol.spgp.h5"
    restart_name = "system_restart.spgr.h5"
    report: dict[str, object] = {}

    with h5py.File(reference_dir / topology_name) as lhs, h5py.File(candidate_dir / topology_name) as rhs:
        topology_paths = _dataset_paths(lhs)
        if topology_paths != _dataset_paths(rhs):
            raise AssertionError("topology dataset schema differs")
        for path in topology_paths:
            if lhs[path].dtype != rhs[path].dtype:
                raise AssertionError(
                    f"{path} dtype differs: {lhs[path].dtype} and {rhs[path].dtype}"
                )
        for path in (
            "/atoms/residue_index",
            "/residues/atom_offset",
            "/parameters/xponge/atoms/name",
            "/parameters/xponge/atoms/type_name",
            "/parameters/xponge/residues/name",
        ):
            _assert_exact(lhs[path][:], rhs[path][:], path)
        report["mass_max_abs"] = _assert_close(lhs["/atoms/mass"][:], rhs["/atoms/mass"][:], "mass")
        report["charge_max_abs"] = _assert_close(lhs["/atoms/charge"][:], rhs["/atoms/charge"][:], "charge")
        report["bond_max_abs"] = _compare_records(lhs, rhs, "/forcefield/bond", ("k", "r0"))
        report["angle_max_abs"] = _compare_records(lhs, rhs, "/forcefield/angle", ("k", "theta0"))
        report["dihedral_max_abs"] = _compare_dihedrals(lhs, rhs)
        report["nb14_max_abs"] = _compare_records(lhs, rhs, "/forcefield/nb14", ("params",))
        report["soft_bond_max_abs"] = _compare_records(
            lhs, rhs, "/forcefield/bond_soft", ("k", "r0", "from_a_or_b")
        )
        listed_root = "/forcefield/custom_force/listed"
        listed_data = f"{listed_root}/data/Ryckaert_Bellemans"
        for path in (
            f"{listed_root}/name",
            f"{listed_root}/potential",
            f"{listed_root}/parameters/text",
            f"{listed_root}/parameters/type",
            f"{listed_root}/parameters/name",
            f"{listed_root}/parameters/offset",
            f"{listed_root}/parameters/is_atom",
            f"{listed_root}/parameters/is_int",
            f"{listed_root}/connected_atoms",
            f"{listed_root}/constrain_distance",
            f"{listed_root}/count",
            f"{listed_data}/name",
            f"{listed_data}/item_count",
            f"{listed_data}/parameter/name",
            f"{listed_data}/parameter/type",
            f"{listed_data}/parameter/is_int",
            f"{listed_data}/parameter/int_value",
        ):
            _assert_exact(lhs[path][()], rhs[path][()], path)
        report["ryckaert_bellemans_max_abs"] = _assert_close(
            lhs[f"{listed_data}/parameter/value"][:],
            rhs[f"{listed_data}/parameter/value"][:],
            "Ryckaert-Bellemans values",
        )
        _assert_close(
            lhs[f"{listed_data}/parameter/float_value"][:],
            rhs[f"{listed_data}/parameter/float_value"][:],
            "Ryckaert-Bellemans typed float values",
        )
        report["lj_max_abs"] = _assert_close(_expanded_lj(lhs), _expanded_lj(rhs), "expanded LJ")
        if _exclusion_pairs(lhs) != _exclusion_pairs(rhs):
            raise AssertionError("exclusion pair sets differ")
        report["atom_count"] = int(lhs["/topology/atom_count"][()])
        report["exclusion_pair_count"] = len(_exclusion_pairs(lhs))

    with h5py.File(reference_dir / restart_name) as lhs, h5py.File(candidate_dir / restart_name) as rhs:
        report["position_max_abs"] = _assert_close(
            lhs["/particles/all/position/value"][:],
            rhs["/particles/all/position/value"][:],
            "restart position",
        )
        report["box_max_abs"] = _assert_close(
            lhs["/particles/all/box/edges/value"][:],
            rhs["/particles/all/box/edges/value"][:],
            "restart box",
        )

    for directory in (reference_dir, candidate_dir):
        with h5py.File(directory / topology_name) as topology, h5py.File(directory / protocol_name) as protocol:
            topology_hash = topology["/topology/topology_hash"].asstr()[()]
            compatible_hash = protocol["/protocol/topology_compatibility/topology_hash"].asstr()[()]
            if not topology_hash or topology_hash != compatible_hash:
                raise AssertionError(f"inconsistent topology compatibility hash in {directory}")
            if topology_hash != _content_hash(topology, "topology.spgt.h5"):
                raise AssertionError(f"invalid topology dataset SHA-256 in {directory}")
            if topology["/topology/atom_order_hash"].asstr()[()] != _content_hash(
                topology, "topology.spgt.h5", ("/atoms/", "/residues/")
            ):
                raise AssertionError(f"invalid atom-order dataset SHA-256 in {directory}")
            if topology["/topology/forcefield_hash"].asstr()[()] != _content_hash(
                topology, "topology.spgt.h5", ("/forcefield/", "/manybody/", "/qc/")
            ):
                raise AssertionError(f"invalid force-field dataset SHA-256 in {directory}")
            if protocol["/identity/content_hash"].asstr()[()] != _content_hash(
                protocol, "protocol.spgp.h5"
            ):
                raise AssertionError(f"invalid protocol dataset SHA-256 in {directory}")
            _assert_exact(protocol["/protocol/cv_count"][()], 0, "protocol CV count")
            _assert_exact(protocol["/protocol/restraint_count"][()], 0, "protocol restraint count")

    report["result"] = "PASS"
    report["hash_semantics"] = "xponge-dataset-sha256"
    report["representational_differences"] = [
        "interaction row order",
        "equivalent LJ type compression",
        "content hashes reflect those dataset representation differences",
    ]
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xponge-root", type=Path, default=Path(__file__).resolve().parents[2] / "XPONGE")
    parser.add_argument("--work-dir", type=Path, help="keep generated bundles in this directory")
    args = parser.parse_args()

    xponge_root = args.xponge_root.resolve()
    if not (xponge_root / "Xponge" / "__init__.py").is_file():
        parser.error(f"Xponge source tree not found: {xponge_root}")

    temporary = None
    if args.work_dir is None:
        temporary = tempfile.mkdtemp(prefix="xponge-bundle-ab-")
        work_dir = Path(temporary)
    else:
        work_dir = args.work_dir.resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
    reference_dir = work_dir / "xponge"
    candidate_dir = work_dir / "xpongecpp"
    reference_dir.mkdir(exist_ok=True)
    candidate_dir.mkdir(exist_ok=True)

    try:
        _run_builder("xponge", reference_dir, xponge_root)
        _run_builder("xpongecpp", candidate_dir, xponge_root)
        print(json.dumps(compare_bundles(reference_dir, candidate_dir), indent=2, sort_keys=True))
        return 0
    finally:
        if temporary is not None:
            shutil.rmtree(temporary)


if __name__ == "__main__":
    raise SystemExit(main())
