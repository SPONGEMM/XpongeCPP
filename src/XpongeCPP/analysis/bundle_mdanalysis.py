"""Expose validated SPONGE input bundles as MDAnalysis universes."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..io_bundle import BundleCase, BundleReader


def _require_mdanalysis():
    try:
        import MDAnalysis as mda
    except ImportError as exc:  # pragma: no cover - required runtime dependency
        raise ModuleNotFoundError(
            "'MDAnalysis' is required to load SPONGE bundles"
        ) from exc
    return mda


def _decode_text_array(values) -> list[str]:
    decoded = []
    for value in np.asarray(values).reshape(-1):
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        elif isinstance(value, np.bytes_):
            decoded.append(bytes(value).decode("utf-8"))
        else:
            decoded.append(str(value))
    return decoded


def _dimensions_from_edges(edges) -> np.ndarray:
    cell = np.asarray(edges, dtype=np.float64).reshape(3, 3)
    lengths = np.linalg.norm(cell, axis=1)
    if np.any(lengths <= 0):
        raise ValueError("bundle box edges must have positive lengths")

    def angle(lhs, rhs, lhs_length, rhs_length):
        cosine = np.dot(lhs, rhs) / (lhs_length * rhs_length)
        return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))

    alpha = angle(cell[1], cell[2], lengths[1], lengths[2])
    beta = angle(cell[0], cell[2], lengths[0], lengths[2])
    gamma = angle(cell[0], cell[1], lengths[0], lengths[1])
    return np.asarray(
        [lengths[0], lengths[1], lengths[2], alpha, beta, gamma],
        dtype=np.float32,
    )


def _case_from_paths(
    topology: str | Path,
    trajectory: str | Path | None,
    protocol: str | Path | None,
    *,
    particle_stream: str,
) -> BundleCase:
    topology_path = Path(topology).resolve()
    root = topology_path.parent
    protocol_path = Path(protocol).resolve() if protocol is not None else None
    restart_path = None
    trajectory_path = None
    if trajectory is not None:
        trajectory_candidate = Path(trajectory).resolve()
        if trajectory_candidate.name.endswith(".spgr.h5"):
            restart_path = trajectory_candidate
        else:
            trajectory_path = trajectory_candidate

    stem_suffix = "_topology.spgt.h5"
    if topology_path.name.endswith(stem_suffix):
        prefix = topology_path.name[: -len(stem_suffix)]
        if protocol_path is None:
            candidate = root / f"{prefix}_protocol.spgp.h5"
            if candidate.is_file():
                protocol_path = candidate
        if restart_path is None and trajectory_path is None:
            candidate = root / f"{prefix}_restart.spgr.h5"
            if candidate.is_file():
                restart_path = candidate

    return BundleCase(
        root=root,
        mdin_path=None,
        mdin_text="",
        commands={},
        topology_path=topology_path,
        protocol_path=protocol_path,
        restart_path=restart_path,
        trajectory_path=trajectory_path,
        restart_load="structural",
        particle_stream=particle_stream,
    )


def load_bundle_universe(
    topology: BundleCase | str | Path,
    trajectory: str | Path | None = None,
    *,
    protocol: str | Path | None = None,
    particle_stream: str = "all",
    frame: int = -1,
    strict: bool = True,
):
    """Load a validated topology plus restart/H5MD frame into MDAnalysis.

    ``topology`` may be a :class:`BundleCase` or a topology HDF5 path. When a
    canonical ``*_topology.spgt.h5`` path is supplied, sibling protocol and
    restart files are discovered automatically. The current implementation
    materializes one selected frame in memory; streaming output trajectories
    remain a later compatibility milestone.
    """

    case = (
        topology
        if isinstance(topology, BundleCase)
        else _case_from_paths(
            topology,
            trajectory,
            protocol,
            particle_stream=particle_stream,
        )
    )
    mda = _require_mdanalysis()
    with BundleReader(case, strict=strict) as reader:
        atom_count = int(
            reader.read_scalar("topology.spgt.h5", "/topology/atom_count")
        )
        residue_index = np.asarray(
            reader.read("topology.spgt.h5", "/atoms/residue_index"),
            dtype=np.int32,
        )
        if residue_index.shape != (atom_count,):
            raise ValueError(
                f"bundle residue_index has shape {residue_index.shape}, "
                f"expected ({atom_count},)"
            )
        residue_count = (
            int(residue_index.max()) + 1 if residue_index.size else 0
        )
        universe = mda.Universe.empty(
            atom_count,
            n_residues=residue_count,
            n_segments=1,
            atom_resindex=residue_index,
            residue_segindex=np.zeros(residue_count, dtype=np.int32),
            trajectory=True,
        )

        masses = np.asarray(
            reader.read("topology.spgt.h5", "/atoms/mass"), dtype=np.float64
        )
        charges = (
            np.asarray(
                reader.read("topology.spgt.h5", "/atoms/charge"),
                dtype=np.float64,
            )
            / 18.2223
        )
        universe.add_TopologyAttr("masses", masses)
        universe.add_TopologyAttr("charges", charges)
        if reader.contains(
            "topology.spgt.h5", "/parameters/xponge/atoms/name"
        ):
            universe.add_TopologyAttr(
                "names",
                _decode_text_array(
                    reader.read(
                        "topology.spgt.h5",
                        "/parameters/xponge/atoms/name",
                    )
                ),
            )
        if reader.contains(
            "topology.spgt.h5", "/parameters/xponge/atoms/type_name"
        ):
            universe.add_TopologyAttr(
                "types",
                _decode_text_array(
                    reader.read(
                        "topology.spgt.h5",
                        "/parameters/xponge/atoms/type_name",
                    )
                ),
            )
        if reader.contains(
            "topology.spgt.h5", "/parameters/xponge/residues/name"
        ):
            universe.add_TopologyAttr(
                "resnames",
                _decode_text_array(
                    reader.read(
                        "topology.spgt.h5",
                        "/parameters/xponge/residues/name",
                    )
                ),
            )
        universe.add_TopologyAttr(
            "ids", np.arange(1, atom_count + 1, dtype=np.int32)
        )
        universe.add_TopologyAttr(
            "resids", np.arange(1, residue_count + 1, dtype=np.int32)
        )
        universe.add_TopologyAttr(
            "resnums", np.arange(1, residue_count + 1, dtype=np.int32)
        )
        universe.add_TopologyAttr("segids", ["SYSTEM"])

        for attr_name, dataset_path in (
            ("bonds", "/forcefield/bond/atoms"),
            ("angles", "/forcefield/angle/atoms"),
            ("dihedrals", "/forcefield/dihedral/atoms"),
            ("impropers", "/forcefield/improper/atoms"),
        ):
            if reader.contains("topology.spgt.h5", dataset_path):
                values = np.asarray(
                    reader.read("topology.spgt.h5", dataset_path),
                    dtype=np.int32,
                )
                if values.size:
                    universe.add_TopologyAttr(attr_name, values)

        bundle_file = (
            "restart.spgr.h5"
            if reader.has_bundle_file("restart.spgr.h5")
            else "trajectory.spg.h5md"
        )
        stream = "all" if bundle_file == "restart.spgr.h5" else case.particle_stream
        position_path = f"/particles/{stream}/position/value"
        positions = np.asarray(
            reader.read(bundle_file, position_path), dtype=np.float32
        )
        selected_frame = frame if frame >= 0 else positions.shape[0] + frame
        if selected_frame < 0 or selected_frame >= positions.shape[0]:
            raise IndexError(
                f"bundle frame {frame} is outside 0..{positions.shape[0] - 1}"
            )
        universe.atoms.positions = positions[selected_frame]

        box_path = f"/particles/{stream}/box/edges/value"
        if reader.contains(bundle_file, box_path):
            edges = np.asarray(reader.read(bundle_file, box_path))
            universe.dimensions = _dimensions_from_edges(edges[selected_frame])
        return universe


__all__ = ["load_bundle_universe"]
