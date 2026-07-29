from pathlib import Path
import runpy
import subprocess

import numpy as np
import pytest

import XpongeCPP as Xponge
from XpongeCPP.io_bundle import (
    BundleCapabilityError,
    BundleConflictError,
    BundleReader,
    BundleValidationError,
    canonical_dataset_hash,
    convert_bundle_to_legacy,
    convert_legacy_to_bundle,
    scan_bundle_case,
)
from io_bundle_fixtures import find_sponge_executable


def _peptide():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    molecule = Xponge.get_peptide_from_sequence("AA")
    molecule.set_box_padding(12.0)
    return molecule


def _legacy_case(tmp_path: Path):
    raw_dir = tmp_path / "raw"
    molecule = _peptide()
    Xponge.save_sponge_input(molecule, "system", raw_dir, format="raw")
    keys = (
        "residue",
        "resname",
        "atom_name",
        "atom_type_name",
        "mass",
        "charge",
        "coordinate",
        "LJ",
        "bond",
        "angle",
        "dihedral",
        "exclude",
        "nb14",
    )
    mdin = "\n".join(
        [
            'mode = "minimization"',
            "step_limit = 0",
            "cutoff = 8.0",
            "print_zeroth_frame = 1",
            "write_mdout_interval = 1",
            "write_information_interval = 1",
            *[
                f'{key}_in_file = "system_{key}.txt"'
                for key in keys
            ],
        ]
    )
    (raw_dir / "mdin.spg.toml").write_text(mdin + "\n")
    return raw_dir


def _numeric_tokens(path: Path):
    return np.asarray(
        [float(token) for token in path.read_text().split()], dtype=np.float64
    )


def _lj_atom_pair_tables(path: Path):
    values = _numeric_tokens(path)
    atom_count, type_count = (int(values[0]), int(values[1]))
    pair_count = type_count * (type_count + 1) // 2
    pair_a = values[2 : 2 + pair_count]
    pair_b = values[2 + pair_count : 2 + pair_count * 2]
    atom_types = values[2 + pair_count * 2 :].astype(np.int64)
    assert atom_types.shape == (atom_count,)

    def expand(values):
        table = np.zeros((type_count, type_count))
        cursor = 0
        for row in range(type_count):
            width = row + 1
            table[row, :width] = values[cursor : cursor + width]
            table[:width, row] = values[cursor : cursor + width]
            cursor += width
        return table[np.ix_(atom_types, atom_types)]

    return expand(pair_a), expand(pair_b)


def _read_hash_inputs(handle):
    excluded = {
        "/schema/name",
        "/schema/version",
        "/parameters/sponge/schema/name",
        "/parameters/sponge/schema/version",
        "/topology/atom_count",
        "/topology/atom_order_hash",
        "/topology/topology_hash",
        "/topology/forcefield_hash",
        "/identity/uuid",
    }
    values = {}

    def collect(name, item):
        path = "/" + name
        if not hasattr(item, "dtype") or path in excluded:
            return
        if item.dtype.kind in {"O", "S", "U"}:
            values[path] = np.asarray(item.asstr()[...], dtype=object)
        else:
            values[path] = np.asarray(item[...])

    handle.visititems(collect)
    return values


def _read_mdout_first_frame(path: Path):
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise AssertionError(f"SPONGE did not write a frame to {path}")
    headers = lines[0].split()
    values = lines[1].split()
    return {
        key: float(value)
        for key, value in zip(headers, values)
        if key not in {"step"}
    }


def test_convert_legacy_to_bundle_round_trips_supported_raw_case(tmp_path):
    h5py = pytest.importorskip("h5py")
    raw_dir = _legacy_case(tmp_path)
    converted_root = tmp_path / "converted"

    manifest = convert_legacy_to_bundle(raw_dir, converted_root)
    bundle_dir = converted_root / "bundle"
    case = scan_bundle_case(bundle_dir)
    with BundleReader(case) as reader:
        atom_count = int(
            reader.read_scalar("topology.spgt.h5", "/topology/atom_count")
        )
        assert atom_count == _peptide().atom_count
    with h5py.File(bundle_dir / "topology.spgt.h5", "r") as handle:
        hash_inputs = _read_hash_inputs(handle)
        expected = canonical_dataset_hash(
            "topology.spgt.h5", hash_inputs
        )
        assert handle["/topology/topology_hash"].asstr()[()] == expected
        parity_hash = runpy.run_path(
            Path(__file__).parents[1]
            / "scripts"
            / "validate_bundle_parity.py"
        )["_content_hash"]
        assert handle["/topology/topology_hash"].asstr()[()] == parity_hash(
            handle, "topology.spgt.h5"
        )

    roundtrip = tmp_path / "roundtrip"
    convert_bundle_to_legacy(bundle_dir, roundtrip)
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
        assert _numeric_tokens(roundtrip / f"input_{key}.txt") == pytest.approx(
            _numeric_tokens(raw_dir / f"system_{key}.txt"),
            abs=2e-5,
            rel=2e-6,
        )
    expected_lj = _lj_atom_pair_tables(raw_dir / "system_LJ.txt")
    actual_lj = _lj_atom_pair_tables(roundtrip / "input_LJ.txt")
    assert actual_lj[0] == pytest.approx(expected_lj[0])
    assert actual_lj[1] == pytest.approx(expected_lj[1])
    for key in ("resname", "atom_name", "atom_type_name"):
        assert (roundtrip / f"input_{key}.txt").read_text() == (
            raw_dir / f"system_{key}.txt"
        ).read_text()
    assert manifest.bundled_mdin == str(
        bundle_dir / "mdin.bundled.spg.toml"
    )


def test_legacy_to_bundle_dry_run_has_no_filesystem_side_effect(tmp_path):
    _h5py = pytest.importorskip("h5py")
    raw_dir = _legacy_case(tmp_path)
    output_dir = tmp_path / "dry"

    manifest = convert_legacy_to_bundle(
        raw_dir, output_dir, dry_run=True
    )

    assert not output_dir.exists()
    assert manifest.schema == "xponge.legacy_to_bundle.manifest"


def test_legacy_to_bundle_rejects_existing_output_before_writing(tmp_path):
    _h5py = pytest.importorskip("h5py")
    raw_dir = _legacy_case(tmp_path)
    output_dir = tmp_path / "existing"
    output_dir.mkdir()
    sentinel = output_dir / "keep.txt"
    sentinel.write_text("keep\n")

    with pytest.raises(BundleConflictError, match="already exists"):
        convert_legacy_to_bundle(raw_dir, output_dir)

    assert sentinel.read_text() == "keep\n"


def test_legacy_to_bundle_rejects_malformed_counted_input(tmp_path):
    _h5py = pytest.importorskip("h5py")
    raw_dir = _legacy_case(tmp_path)
    (raw_dir / "system_mass.txt").write_text("2\n12.0\n")

    with pytest.raises(BundleValidationError, match="declares 2 mass"):
        convert_legacy_to_bundle(raw_dir, tmp_path / "converted")


def test_legacy_to_bundle_rejects_unknown_present_input(tmp_path):
    _h5py = pytest.importorskip("h5py")
    raw_dir = _legacy_case(tmp_path)
    (raw_dir / "custom.txt").write_text("payload\n")
    with (raw_dir / "mdin.spg.toml").open("a") as handle:
        handle.write('custom_in_file = "custom.txt"\n')

    with pytest.raises(BundleCapabilityError, match="custom_in_file"):
        convert_legacy_to_bundle(raw_dir, tmp_path / "converted")


def test_legacy_converter_is_available_from_legacy_package():
    from Xponge.io_bundle.converter import (
        LegacyToBundleConverter as LegacyImport,
    )
    from XpongeCPP.io_bundle import LegacyToBundleConverter

    assert LegacyImport is LegacyToBundleConverter


def test_round_tripped_legacy_inputs_match_sponge_zeroth_frame(tmp_path):
    executable = find_sponge_executable()
    if executable is None:
        pytest.skip(
            "SPONGE was not configured, found on PATH, or built in the "
            "sibling SPONGE checkout"
        )
    raw_dir = _legacy_case(tmp_path)
    converted_root = tmp_path / "converted"
    convert_legacy_to_bundle(raw_dir, converted_root)
    roundtrip_dir = tmp_path / "roundtrip"
    convert_bundle_to_legacy(converted_root / "bundle", roundtrip_dir)

    runs = (
        (raw_dir, "mdin.spg.toml"),
        (roundtrip_dir, "mdin.legacy.spg.toml"),
    )
    frames = []
    for cwd, mdin in runs:
        completed = subprocess.run(
            [str(executable), "-mdin", mdin],
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        assert completed.returncode == 0, (
            f"SPONGE failed in {cwd}\nstdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
        frames.append(_read_mdout_first_frame(cwd / "mdout.txt"))

    assert frames[0].keys() == frames[1].keys()
    for key in frames[0]:
        assert frames[1][key] == pytest.approx(
            frames[0][key], abs=1e-5, rel=1e-7
        ), key
