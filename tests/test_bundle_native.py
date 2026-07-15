import pytest

import XpongeCPP as Xponge


def _peptide():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.get_peptide_from_sequence("AA")


def test_native_bundle_saver_writes_typed_hdf5(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()

    returned = Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    paths = {
        "topology": tmp_path / "system_topology.spgt.h5",
        "protocol": tmp_path / "system_protocol.spgp.h5",
        "restart": tmp_path / "system_restart.spgr.h5",
    }

    assert returned is molecule
    with h5py.File(paths["topology"], "r") as handle:
        assert handle["/schema/name"].asstr()[()] == "sponge.topology.h5"
        assert handle["/atoms/mass"].shape == (molecule.atom_count,)
        assert handle["/atoms/charge"].shape == (molecule.atom_count,)
        assert handle["/atoms/residue_index"].shape == (molecule.atom_count,)
        assert handle["/forcefield/bond/atoms"].shape[1] == 2
        assert handle["/forcefield/angle/atoms"].shape[1] == 3
        assert handle["/forcefield/dihedral/atoms"].shape[1] == 4
        assert "/forcefield/improper" not in handle
        assert handle["/forcefield/nb14/params"].shape[1] == 2
        assert int(handle["/topology/atom_count"][()]) == molecule.atom_count
        assert handle["/atoms/charge"].attrs["unit"] == "Amber"
        offsets = handle["/topology/exclusions/offset"][:]
        exclusions = handle["/topology/exclusions/list"][:]
        for atom_index in range(molecule.atom_count):
            row = exclusions[offsets[atom_index] : offsets[atom_index + 1]]
            assert all(atom_index < excluded_atom for excluded_atom in row)
    with h5py.File(paths["protocol"], "r") as handle:
        assert handle["/schema/name"].asstr()[()] == "sponge.protocol.h5"
    with h5py.File(paths["restart"], "r") as handle:
        assert handle["/particles/all/position/value"].shape == (1, molecule.atom_count, 3)
        assert handle["/particles/all/box/edges/value"].shape == (1, 3, 3)
        assert tuple(handle["/h5md"].attrs["version"]) == (1, 1)
        assert handle["/particles/all/position/value"].attrs["unit"] == "Angstrom"
        assert handle["/particles/all/position/step"].id == handle["/particles/all/step"].id


def test_native_bundle_saver_is_available_from_legacy_package(tmp_path):
    import Xponge

    molecule = _peptide()
    assert Xponge.save_sponge_input_bundle(molecule, "legacy", tmp_path) is molecule
    assert (tmp_path / "legacy_topology.spgt.h5").is_file()


def test_native_bundle_saver_accepts_residue_type(tmp_path):
    _peptide()
    residue_type = Xponge.ResidueType.get_type("ALA")

    molecule = Xponge.save_sponge_input_bundle(residue_type, "residue_type", tmp_path)

    assert molecule.residue_count == 1
    assert (tmp_path / "residue_type_topology.spgt.h5").is_file()


def test_native_bundle_saver_materializes_extended_topology(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()
    molecule.add_virtual_atom2(0, 1, 2, 3, 0.25, 0.75)
    molecule.add_nb14_extra(4, 1, 1.25, 2.5, 0.75)
    molecule.add_urey_bradley(2, 1, 0, 1.1, 2.2, 3.3, 4.4)
    cmap_type = molecule.add_cmap_type(2, [0.1, 0.2, 0.3, 0.4])
    molecule.add_cmap(4, 3, 2, 1, 0, cmap_type)
    molecule.set_gb_radius("bondi_radii")
    molecule.enable_subsys_division()

    Xponge.save_sponge_input_bundle(molecule, "extended", tmp_path)

    with h5py.File(tmp_path / "extended_topology.spgt.h5", "r") as handle:
        assert handle["/forcefield/virtual_atom/from"].shape == (3,)
        assert handle["/forcefield/nb14_extra/params"].shape == (1, 3)
        assert handle["/forcefield/urey_bradley/atoms"].shape == (1, 3)
        assert handle["/forcefield/cmap/atoms"].shape == (1, 5)
        assert handle["/forcefield/cmap/grid_value"].shape == (4,)
        assert handle["/forcefield/gb/params"].shape == (molecule.atom_count, 2)
        assert handle["/forcefield/subsys_division"].shape == (molecule.atom_count,)


def test_native_bundle_saver_uses_sponge_improper_atom_order(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()
    molecule.add_improper_dihedral(1, 2, 3, 4, 5.5, 180.0)

    Xponge.save_sponge_input_bundle(molecule, "improper", tmp_path)

    with h5py.File(tmp_path / "improper_topology.spgt.h5", "r") as handle:
        assert handle["/forcefield/improper/atoms"][-1].tolist() == [3, 1, 2, 4]


def test_native_bundle_saver_materializes_lj_soft_core(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()
    molecule.atoms[0].lj_type_b = molecule.atoms[1].type
    molecule.enable_lj_soft_core()

    Xponge.save_sponge_input_bundle(molecule, "soft", tmp_path)

    with h5py.File(tmp_path / "soft_topology.spgt.h5", "r") as handle:
        assert "/forcefield/lj" not in handle
        assert handle["/forcefield/lj_soft_core/atom_type_A"].shape == (molecule.atom_count,)
        assert handle["/forcefield/lj_soft_core/atom_type_B"].shape == (molecule.atom_count,)
        assert handle["/forcefield/lj_soft_core/pair_AA"].ndim == 1
        assert handle["/forcefield/lj_soft_core/pair_BB"].ndim == 1


def test_native_bundle_saver_rejects_escaping_prefix(tmp_path):
    with pytest.raises(ValueError, match="escapes output directory"):
        Xponge.save_sponge_input_bundle(_peptide(), "../escape", tmp_path)
