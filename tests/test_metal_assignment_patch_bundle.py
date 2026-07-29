"""XpongeCPP-native integration gates for the origin 1.7b8 patch API."""

from __future__ import annotations

import numpy as np

import XpongeCPP as Xponge
from XpongeCPP.metal_assignment import apply

from origin_test_metal_assignment_apply import _local_patch, _ordinary_molecule


def test_local_patch_exports_the_same_atom_state_to_raw_and_bundle(tmp_path):
    h5py = __import__("h5py")
    patch = _local_patch()
    molecule, mapping = _ordinary_molecule(patch)
    applied = apply(molecule, patch, mapping).molecule
    applied.set_box_padding(20.0)

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    bundle_dir = tmp_path / "bundle"
    bundle_dir.mkdir()
    Xponge.save_sponge_input(applied, "system", raw_dir, format="raw")
    Xponge.save_sponge_input(applied, "system", bundle_dir, format="bundle")

    raw_charge_lines = (
        raw_dir / "system_charge.txt"
    ).read_text(encoding="utf-8").splitlines()
    raw_charges = np.asarray(
        [float(value) for value in raw_charge_lines[1:]],
        dtype=float,
    )
    with h5py.File(
        bundle_dir / "system_topology.spgt.h5",
        "r",
    ) as topology:
        bundle_charges = np.asarray(topology["/atoms/charge"][()])
        np.testing.assert_allclose(
            bundle_charges,
            np.asarray([atom.charge for atom in applied.atoms]) * 18.2223,
            rtol=0.0,
            atol=2.0e-6,
        )
        assert topology["/forcefield/bond/atoms"].shape[0] == 1

    np.testing.assert_allclose(
        raw_charges,
        bundle_charges,
        rtol=0.0,
        atol=1.0e-6,
    )


def test_local_patch_inplace_is_a_copy_then_commit():
    patch = _local_patch()
    molecule, mapping = _ordinary_molecule(patch)

    result = apply(molecule, patch, mapping, inplace=True)

    assert result.molecule is molecule
    metal_id = patch.target_metal_atom_ids[0]
    assert molecule.atoms[mapping[metal_id]].charge == (
        patch.parameterization_result.metal_overlay.charges[metal_id]
    )
