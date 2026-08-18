from __future__ import annotations

import numpy as np
import pytest


pytest.importorskip("MDAnalysis")


def test_sponge_input_reader_exposes_origin_format_contract():
    from Xponge.analysis import md_analysis as xmda

    assert xmda.SpongeInputReader.format == "SPONGE_MASS"
    assert xmda.SpongeInputReader._format_hint("system_mass.txt") is True
    assert xmda.SpongeInputReader._format_hint("system_charge.txt") is False


def test_legacy_sponge_dat_reader_writer_and_registry_roundtrip(tmp_path):
    import MDAnalysis as mda
    from Xponge.analysis import md_analysis as xmda

    universe = mda.Universe.empty(2, trajectory=True)
    frames = [
        np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        np.array([[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]]),
    ]
    dimensions = [
        np.array([10.0, 11.0, 12.0, 90.0, 90.0, 90.0]),
        np.array([20.0, 21.0, 22.0, 90.0, 90.0, 90.0]),
    ]
    trajectory_path = tmp_path / "trajectory.dat"

    with xmda.SpongeTrajectoryWriter(str(trajectory_path)) as writer:
        for positions, box in zip(frames, dimensions):
            universe.atoms.positions = positions
            universe.dimensions = box
            writer.write(universe)

    assert xmda.SpongeTrajectoryReader.format == "SPONGE_TRAJ"
    assert xmda.SpongeTrajectoryReader._format_hint(str(trajectory_path))
    assert mda._MULTIFRAME_WRITERS["SPONGE_TRAJ"] is xmda.SpongeTrajectoryWriter
    assert mda._MULTIFRAME_WRITERS["DAT"] is xmda.SpongeTrajectoryWriter

    reader = xmda.SpongeTrajectoryReader(
        str(trajectory_path), n_atoms=2, box=str(tmp_path / "trajectory.box")
    )
    assert reader.n_frames == 2
    np.testing.assert_allclose(reader.ts.positions, frames[0])
    np.testing.assert_allclose(reader.ts.dimensions, dimensions[0])
    np.testing.assert_allclose(reader._read_frame(1).positions, frames[1])
    np.testing.assert_allclose(reader.ts.dimensions, dimensions[1])
    reader.close()

    topology_path = tmp_path / "system_mass.txt"
    topology_path.write_text("2\n12.011\n1.008\n", encoding="utf-8")
    loaded = mda.Universe(
        str(topology_path),
        str(trajectory_path),
        format=xmda.SpongeTrajectoryReader,
        box=str(tmp_path / "trajectory.box"),
    )
    assert loaded.atoms.n_atoms == 2
    np.testing.assert_allclose(loaded.trajectory[1].positions, frames[1])
    np.testing.assert_allclose(loaded.trajectory.ts.dimensions, dimensions[1])

    auto_detected = mda.Universe(
        str(topology_path), str(trajectory_path), box=str(tmp_path / "trajectory.box")
    )
    assert auto_detected.trajectory.__class__ is xmda.SpongeTrajectoryReader
    np.testing.assert_allclose(auto_detected.trajectory[1].positions, frames[1])


def test_legacy_sponge_coordinate_reader_writer_roundtrip(tmp_path):
    import MDAnalysis as mda
    from Xponge.analysis import md_analysis as xmda

    universe = mda.Universe.empty(2, trajectory=True)
    positions = np.array([[1.25, 2.5, 3.75], [4.0, 5.0, 6.0]])
    dimensions = np.array([15.0, 16.0, 17.0, 90.0, 91.0, 92.0])
    universe.atoms.positions = positions
    universe.dimensions = dimensions
    coordinate_path = tmp_path / "system_coordinate.txt"

    with xmda.SpongeCoordinateWriter(str(coordinate_path)) as writer:
        writer.write(universe)

    assert xmda.SpongeCoordinateReader.format == "SPONGE_CRD"
    assert xmda.SpongeCoordinateReader._format_hint(str(coordinate_path))
    assert mda._SINGLEFRAME_WRITERS["SPONGE_CRD"] is xmda.SpongeCoordinateWriter
    assert mda._SINGLEFRAME_WRITERS["TXT"] is xmda.SpongeCoordinateWriter

    reader = xmda.SpongeCoordinateReader(str(coordinate_path), n_atoms=2)
    assert reader.n_frames == 1
    np.testing.assert_allclose(reader.ts.positions, positions)
    np.testing.assert_allclose(reader.ts.dimensions, dimensions)
    reader.close()
