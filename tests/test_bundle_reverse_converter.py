import json
from pathlib import Path

import numpy as np
import pytest

import XpongeCPP as Xponge
from XpongeCPP.io_bundle import (
    BundleCapabilityError,
    BundleConflictError,
    bundle_case_from_prefix,
    convert_bundle_to_legacy,
)


def _peptide():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    molecule = Xponge.get_peptide_from_sequence("AA")
    molecule.set_box_padding(4.0)
    return molecule


def _write_raw_and_bundle(tmp_path: Path):
    raw_dir = tmp_path / "raw"
    bundle_dir = tmp_path / "bundle"
    molecule = _peptide()
    Xponge.save_sponge_input(molecule, "system", raw_dir, format="raw")
    Xponge.save_sponge_input_bundle(molecule, "system", bundle_dir)
    return raw_dir, bundle_dir


def _numeric_tokens(path: Path):
    return np.asarray(
        [float(token) for token in path.read_text().split()], dtype=np.float64
    )


def _lj_atom_pair_tables(path: Path):
    values = _numeric_tokens(path)
    atom_count, type_count = (int(values[0]), int(values[1]))
    triangular_count = type_count * (type_count + 1) // 2
    pair_a = values[2 : 2 + triangular_count]
    pair_b = values[2 + triangular_count : 2 + 2 * triangular_count]
    atom_types = values[2 + 2 * triangular_count :].astype(np.int64)
    assert atom_types.shape == (atom_count,)

    def expand(values):
        table = np.zeros((type_count, type_count), dtype=np.float64)
        cursor = 0
        for row in range(type_count):
            width = row + 1
            table[row, :width] = values[cursor : cursor + width]
            table[:width, row] = values[cursor : cursor + width]
            cursor += width
        return table[np.ix_(atom_types, atom_types)]

    return expand(pair_a), expand(pair_b)


def test_convert_native_bundle_to_legacy_matches_direct_core_export(tmp_path):
    _h5py = pytest.importorskip("h5py")
    raw_dir, bundle_dir = _write_raw_and_bundle(tmp_path)
    output_dir = tmp_path / "converted"

    manifest = convert_bundle_to_legacy(
        bundle_dir,
        output_dir,
        prefix="system",
    )

    numeric_keys = (
        "residue",
        "mass",
        "charge",
        "coordinate",
        "bond",
        "angle",
        "dihedral",
        "exclude",
        "nb14",
    )
    for key in numeric_keys:
        expected = _numeric_tokens(raw_dir / f"system_{key}.txt")
        actual = _numeric_tokens(output_dir / f"system_{key}.txt")
        assert actual.shape == expected.shape, key
        assert actual == pytest.approx(expected, abs=2e-5, rel=2e-6), key
    expected_lj = _lj_atom_pair_tables(raw_dir / "system_LJ.txt")
    actual_lj = _lj_atom_pair_tables(output_dir / "system_LJ.txt")
    assert actual_lj[0] == pytest.approx(
        expected_lj[0], abs=2e-5, rel=2e-6
    )
    assert actual_lj[1] == pytest.approx(
        expected_lj[1], abs=2e-5, rel=2e-6
    )
    for key in ("resname", "atom_name", "atom_type_name"):
        assert (output_dir / f"system_{key}.txt").read_text() == (
            raw_dir / f"system_{key}.txt"
        ).read_text()

    mdin = (output_dir / "mdin.legacy.spg.toml").read_text()
    assert 'coordinate_in_file = "system_coordinate.txt"' in mdin
    assert 'LJ_in_file = "system_LJ.txt"' in mdin
    assert manifest.generated_mdin == str(
        output_dir / "mdin.legacy.spg.toml"
    )
    manifest_data = json.loads(
        (output_dir / "manifest.bundle_to_legacy.json").read_text()
    )
    assert manifest_data["schema"] == "xponge.bundle_to_legacy.manifest"
    assert {entry["key"] for entry in manifest_data["entries"]} >= {
        "coordinate",
        "LJ",
        "nb14",
    }


def test_bundle_to_legacy_dry_run_writes_nothing(tmp_path):
    _h5py = pytest.importorskip("h5py")
    _raw_dir, bundle_dir = _write_raw_and_bundle(tmp_path)
    output_dir = tmp_path / "dry"

    manifest = convert_bundle_to_legacy(
        bundle_dir,
        output_dir,
        prefix="system",
        dry_run=True,
    )

    assert not output_dir.exists()
    assert manifest.generated_mdin == str(
        output_dir / "mdin.legacy.spg.toml"
    )


def test_bundle_to_legacy_checks_all_conflicts_before_writing(tmp_path):
    _h5py = pytest.importorskip("h5py")
    _raw_dir, bundle_dir = _write_raw_and_bundle(tmp_path)
    output_dir = tmp_path / "conflict"
    output_dir.mkdir()
    sentinel = output_dir / "system_mass.txt"
    sentinel.write_text("keep\n")

    with pytest.raises(BundleConflictError, match="already exist"):
        convert_bundle_to_legacy(
            bundle_dir,
            output_dir,
            prefix="system",
        )

    assert sentinel.read_text() == "keep\n"
    assert not (output_dir / "system_charge.txt").exists()


def test_bundle_to_legacy_exports_supported_typed_bond_soft_force(tmp_path):
    _h5py = pytest.importorskip("h5py")
    bundle_dir = tmp_path / "bundle"
    molecule = _peptide()
    molecule.add_bond_soft(1, 0, 12.5, 1.25, 1)
    Xponge.save_sponge_input_bundle(molecule, "soft", bundle_dir)

    converted = tmp_path / "converted"
    convert_bundle_to_legacy(
        bundle_dir,
        converted,
        prefix="soft",
    )
    values = (converted / "soft_bond_soft.txt").read_text().split()
    assert values == ["1", "0", "1", "12.5", "1.25", "1"]


def test_bundle_to_legacy_class_and_legacy_import_path(tmp_path):
    _h5py = pytest.importorskip("h5py")
    bundle_dir = tmp_path / "bundle"
    molecule = _peptide()
    Xponge.save_sponge_input_bundle(molecule, "system", bundle_dir)
    case = bundle_case_from_prefix(bundle_dir, "system")

    from Xponge.io_bundle.reverse_converter import (
        BundleToLegacyConverter as LegacyConverter,
    )
    from XpongeCPP.io_bundle import BundleToLegacyConverter

    assert LegacyConverter is BundleToLegacyConverter
    manifest = BundleToLegacyConverter(
        case, tmp_path / "converted", prefix="system"
    ).convert(dry_run=True)
    assert manifest.schema_version == 1
