"""Direct XpongeCPP molecule-to-bundle saver with native protocol support."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np

from .bundle_builder import (
    BundleBuilder,
    BundleMetadata,
    BundlePaths,
    _RESTART_PARTICLE_STATE_PATHS,
)
from .errors import BundlePathError
from .protocol import SpongeProtocol, add_protocol_to_bundle


def save_sponge_input_bundle(
    molecule,
    prefix=None,
    dirname=".",
    *,
    protocol: SpongeProtocol | None = None,
    source_atom_ids=None,
    return_mapping=False,
):
    """Write native C++ topology/restart data plus an origin-compatible protocol."""

    from .. import _core, get_template_molecule, has_template

    target = molecule
    if isinstance(molecule, _core.Residue):
        target = _core.Molecule(molecule.name)
        target.add_residue(molecule)
        if target.atom_count != molecule.atom_count:
            raise ValueError("residue materialization changed atom coverage")
        for source_atom, target_atom in zip(molecule.atoms, target.atoms):
            for field in (
                "x",
                "y",
                "z",
                "charge",
                "mass",
                "type",
                "bad_coordinate",
            ):
                if hasattr(source_atom, field) and hasattr(target_atom, field):
                    setattr(target_atom, field, getattr(source_atom, field))
        coordinates = np.asarray(
            [[atom.x, atom.y, atom.z] for atom in target.atoms], dtype=float
        )
        lengths = np.maximum(np.abs(coordinates).max(axis=0) + 5.0, 10.0)
        target.set_periodic_box([0.0, 0.0, 0.0], lengths.tolist())
    elif isinstance(molecule, _core.ResidueType):
        target = _core.molecule_from_residuetype(molecule)
    elif not isinstance(molecule, _core.Molecule):
        if hasattr(molecule, "name") and has_template(molecule.name):
            target = get_template_molecule(molecule.name)
        else:
            raise TypeError(
                "save_sponge_input_bundle expects a Molecule, Residue, "
                "ResidueType, or registered template-like object"
            )

    output_root = Path(dirname).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    normalized_prefix = str(prefix or getattr(molecule, "name", "system"))
    final_paths = _prefixed_bundle_paths(output_root, normalized_prefix)

    with tempfile.TemporaryDirectory(
        prefix=".xpongecpp-bundle-", dir=output_root
    ) as tempdir:
        staging_root = Path(tempdir)
        staged_prefix = "staged"
        prepared = _core.save_sponge_input_bundle(
            target, staged_prefix, str(staging_root), protocol=None
        )
        staged_paths = _prefixed_bundle_paths(staging_root, staged_prefix)
        _apply_protocol(staged_paths, protocol)
        for source, destination in (
            (staged_paths.topology, final_paths.topology),
            (staged_paths.protocol, final_paths.protocol),
            (staged_paths.restart, final_paths.restart),
        ):
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, destination)
    del prepared
    if not return_mapping:
        return target
    if source_atom_ids is None:
        raise ValueError("return_mapping=True requires source_atom_ids")
    atoms = list(target.atoms)
    if isinstance(source_atom_ids, dict):
        by_index = {
            int(atom.index): str(value) for atom, value in source_atom_ids.items()
        }
        if set(by_index) != {int(atom.index) for atom in atoms}:
            raise ValueError(
                "source_atom_ids mapping must cover every input Atom exactly once"
            )
        values = tuple(by_index[int(atom.index)] for atom in atoms)
    else:
        values = tuple(str(value) for value in source_atom_ids)
        if len(values) != len(atoms):
            raise ValueError(
                "source_atom_ids must contain one ID for every input Atom"
            )
    if len(set(values)) != len(values):
        raise ValueError("source_atom_ids must be unique")
    return target, tuple(
        {"simulation_index": index, "source_atom_id": source_id}
        for index, source_id in enumerate(values)
    )


def _apply_protocol(paths: BundlePaths, protocol: SpongeProtocol | None) -> None:
    import h5py

    with h5py.File(paths.topology, "r") as topology:
        atom_count = int(np.asarray(topology["/topology/atom_count"][()]))
        topology_hash = _text(topology["/topology/topology_hash"][()])
        atom_order_hash = _text(topology["/topology/atom_order_hash"][()])
        forcefield_hash = _text(topology["/topology/forcefield_hash"][()])
        identity_uuid = _text(topology["/identity/uuid"][()])

    builder = BundleBuilder(paths, identity_uuid=identity_uuid)
    with h5py.File(paths.restart, "r") as restart:
        for dataset_path in _RESTART_PARTICLE_STATE_PATHS:
            if dataset_path in restart:
                builder._dataset_hashes[
                    f"restart.spgr.h5:{dataset_path}"
                ] = np.asarray(restart[dataset_path][...])

    summary = add_protocol_to_bundle(protocol, builder, atom_count=atom_count)
    metadata = BundleMetadata(
        atom_count=atom_count,
        topology_hash=topology_hash,
        atom_order_hash=atom_order_hash,
        forcefield_hash=forcefield_hash,
        protocol_hash=builder.content_hash(bundle_file="protocol.spgp.h5"),
        cv_count=summary.cv_count,
        restraint_count=summary.restraint_count,
        enhanced_methods=summary.enhanced_methods,
        creator_version="xpongecpp-native-protocol",
    )
    builder.finalize({"protocol.spgp.h5", "restart.spgr.h5"}, metadata)


def _prefixed_bundle_paths(root: Path, prefix: str) -> BundlePaths:
    relative = Path(prefix)
    if relative.is_absolute():
        raise BundlePathError(f"bundled input prefix must be relative: {prefix}")
    base = (root / relative).resolve()
    if base == root:
        raise BundlePathError(
            f"bundled input prefix must identify a file prefix below the output directory: {prefix}"
        )
    try:
        base.relative_to(root)
    except ValueError as exc:
        raise BundlePathError(
            f"bundled input prefix escapes output directory: {prefix}"
        ) from exc
    return BundlePaths(
        topology=Path(str(base) + "_topology.spgt.h5"),
        protocol=Path(str(base) + "_protocol.spgp.h5"),
        restart=Path(str(base) + "_restart.spgr.h5"),
        trajectory=None,
    )


def _text(value) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    return str(value)
