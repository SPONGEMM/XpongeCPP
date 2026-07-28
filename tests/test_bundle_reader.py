from pathlib import Path

import pytest

import XpongeCPP as Xponge
from XpongeCPP.io_bundle import (
    BundlePathError,
    BundleReader,
    BundleSchemaError,
    BundleValidationError,
    bundle_case_from_prefix,
    scan_bundle_case,
)


def _write_bundle(tmp_path: Path, prefix: str = "system"):
    pytest.importorskip("h5py")
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    molecule = (
        Xponge.ResidueType.get_type("NALA")
        + Xponge.ResidueType.get_type("ALA")
        + Xponge.ResidueType.get_type("CALA")
    )
    molecule.set_box_padding(4.0)
    Xponge.save_sponge_input_bundle(molecule, prefix, tmp_path)
    return molecule, bundle_case_from_prefix(tmp_path, prefix)


def _replace_text_dataset(path: Path, dataset_path: str, value: str):
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "r+") as handle:
        del handle[dataset_path]
        handle.create_dataset(
            dataset_path, data=value, dtype=h5py.string_dtype(encoding="utf-8")
        )


def test_bundle_reader_validates_and_reads_native_three_file_bundle(tmp_path):
    molecule, case = _write_bundle(tmp_path)

    with BundleReader(case) as reader:
        assert reader.warnings == []
        assert reader.has_bundle_file("topology.spgt.h5")
        assert reader.has_bundle_file("protocol.spgp.h5")
        assert reader.has_bundle_file("restart.spgr.h5")
        assert int(
            reader.read_scalar("topology.spgt.h5", "/topology/atom_count")
        ) == molecule.atom_count
        assert reader.read("topology.spgt.h5", "/atoms/mass").shape == (
            molecule.atom_count,
        )
        assert reader.read(
            "restart.spgr.h5", "/particles/all/position/value"
        ).shape == (1, molecule.atom_count, 3)
        identities = {
            reader.read_text(bundle_file, "/identity/uuid")
            for bundle_file in (
                "topology.spgt.h5",
                "protocol.spgp.h5",
                "restart.spgr.h5",
            )
        }
        assert len(identities) == 1


def test_scan_bundle_case_resolves_relative_mdin_bindings(tmp_path):
    _molecule, direct_case = _write_bundle(tmp_path)
    mdin = tmp_path / "mdin.bundled.spg.toml"
    mdin.write_text(
        "\n".join(
            [
                'input_h5_topology_path = "system_topology.spgt.h5"',
                'input_h5_protocol_path = "system_protocol.spgp.h5"',
                'input_h5_restart_path = "system_restart.spgr.h5"',
                'input_h5_restart_load = "structural"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    scanned = scan_bundle_case(tmp_path)
    assert scanned.topology_path == direct_case.topology_path
    assert scanned.protocol_path == direct_case.protocol_path
    assert scanned.restart_path == direct_case.restart_path
    with BundleReader(scanned):
        pass


def test_bundle_reader_rejects_or_warns_on_schema_version_mismatch(tmp_path):
    _molecule, case = _write_bundle(tmp_path)
    _replace_text_dataset(
        case.protocol_path, "/schema/version", "sponge.input.v1"
    )

    with pytest.raises(BundleSchemaError, match="schema version"):
        with BundleReader(case):
            pass

    with BundleReader(case, strict=False) as reader:
        assert any("schema version" in warning for warning in reader.warnings)


def test_bundle_reader_rejects_lineage_mismatch(tmp_path):
    _molecule, case = _write_bundle(tmp_path)
    _replace_text_dataset(
        case.restart_path, "/run/topology_hash", "sha256:" + "0" * 64
    )

    with pytest.raises(BundleValidationError, match="restart topology hash"):
        with BundleReader(case):
            pass


def test_bundle_reader_rejects_atom_dimension_mismatch(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule, case = _write_bundle(tmp_path)
    with h5py.File(case.topology_path, "r+") as handle:
        del handle["/atoms/mass"]
        handle.create_dataset(
            "/atoms/mass", data=[1.0] * (molecule.atom_count - 1)
        )

    with pytest.raises(BundleValidationError, match="/atoms/mass has shape"):
        with BundleReader(case):
            pass


def test_bundle_reader_rejects_legacy_sidecar_path_escape(tmp_path):
    h5py = pytest.importorskip("h5py")
    _molecule, case = _write_bundle(tmp_path)
    with h5py.File(case.topology_path, "r+") as handle:
        group = handle.require_group("/parameters/sponge/files/legacy_sidecars")
        string_type = h5py.string_dtype(encoding="utf-8")
        group.create_dataset("key", data=["custom"], dtype=string_type)
        group.create_dataset("path", data=["../escape.txt"], dtype=string_type)

    with pytest.raises(BundlePathError, match="escapes bundle root"):
        with BundleReader(case):
            pass
