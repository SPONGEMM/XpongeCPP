import numpy as np
import pytest

import XpongeCPP as Xponge
from XpongeCPP.analysis import load_bundle_universe
from XpongeCPP.io_bundle import (
    BundleReader,
    BundleValidationError,
    bundle_case_from_prefix,
)


def _write_bundle(tmp_path, prefix="universe"):
    pytest.importorskip("MDAnalysis")
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    molecule = (
        Xponge.ResidueType.get_type("NALA")
        + Xponge.ResidueType.get_type("ALA")
        + Xponge.ResidueType.get_type("CALA")
    )
    molecule.set_box_padding(4.0)
    Xponge.save_sponge_input_bundle(molecule, prefix, tmp_path)
    return molecule, bundle_case_from_prefix(tmp_path, prefix)


def test_load_bundle_universe_materializes_topology_and_restart_frame(tmp_path):
    molecule, case = _write_bundle(tmp_path)

    universe = load_bundle_universe(case)

    assert universe.atoms.n_atoms == molecule.atom_count
    assert universe.residues.n_residues == molecule.residue_count
    assert universe.segments.n_segments == 1
    assert list(universe.residues.resnames) == [
        residue.name for residue in molecule.residues
    ]
    assert list(universe.atoms.names) == [atom.name for atom in molecule.atoms]
    assert list(universe.atoms.types) == [atom.type for atom in molecule.atoms]
    assert np.allclose(
        universe.atoms.masses,
        [atom.mass for atom in molecule.atoms],
        atol=1e-6,
    )
    assert np.allclose(
        universe.atoms.charges,
        [atom.charge for atom in molecule.atoms],
        atol=1e-6,
    )
    with BundleReader(case) as reader:
        expected_positions = reader.read(
            "restart.spgr.h5", "/particles/all/position/value"
        )[0]
        expected_bonds = reader.read(
            "topology.spgt.h5", "/forcefield/bond/atoms"
        )
    assert np.allclose(universe.atoms.positions, expected_positions)
    assert len(universe.bonds) == len(expected_bonds)
    assert np.all(universe.dimensions[:3] > 0)
    assert np.allclose(universe.dimensions[3:], [90.0, 90.0, 90.0])


def test_load_bundle_universe_discovers_sibling_artifacts_from_topology(tmp_path):
    molecule, case = _write_bundle(tmp_path, "discovered")

    universe = load_bundle_universe(case.topology_path)

    assert universe.atoms.n_atoms == molecule.atom_count
    assert universe.trajectory.n_frames == 1


def test_load_bundle_universe_preserves_reader_lineage_validation(tmp_path):
    h5py = pytest.importorskip("h5py")
    _molecule, case = _write_bundle(tmp_path, "bad_lineage")
    with h5py.File(case.restart_path, "r+") as handle:
        del handle["/run/atom_order_hash"]
        handle.create_dataset(
            "/run/atom_order_hash",
            data="sha256:" + "0" * 64,
            dtype=h5py.string_dtype(encoding="utf-8"),
        )

    with pytest.raises(BundleValidationError, match="atom-order hash"):
        load_bundle_universe(case)


def test_load_bundle_universe_rejects_out_of_range_frame(tmp_path):
    _molecule, case = _write_bundle(tmp_path, "bad_frame")

    with pytest.raises(IndexError, match="outside"):
        load_bundle_universe(case, frame=1)


def test_legacy_xponge_package_exposes_bundle_reader_and_mdanalysis():
    from Xponge.analysis.bundle_mdanalysis import (
        load_bundle_universe as legacy_load_bundle_universe,
    )
    from Xponge.io_bundle import BundleReader as LegacyBundleReader
    from Xponge.io_bundle.bundle_case import BundleCase as LegacyBundleCase

    assert legacy_load_bundle_universe is load_bundle_universe
    assert LegacyBundleReader is BundleReader
    assert LegacyBundleCase.__module__ == "XpongeCPP.io_bundle.case"
