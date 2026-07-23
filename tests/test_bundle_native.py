import hashlib
import runpy
from pathlib import Path

import numpy as np
import pytest

import XpongeCPP as Xponge


def _peptide():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.get_peptide_from_sequence("AA")


def _content_hash(handle, bundle_file, path_prefixes=()):
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
        if hasattr(item, "dtype") and "/" + name not in excluded
        else None
    )
    digest = hashlib.sha256(bundle_file.encode())
    for path, dataset in sorted(datasets.items()):
        if path_prefixes and not path.startswith(path_prefixes):
            continue
        values = dataset[()]
        array = np.asarray(values)
        dtype = dataset.dtype
        if array.dtype.kind == "b":
            array = array.astype(np.uint8, copy=False)
            dtype = array.dtype
        digest.update(b"\0" + path.encode() + b"\0")
        digest.update(str(dtype).encode("ascii") + b"\0")
        digest.update(repr(array.shape).encode("ascii") + b"\0")
        if array.dtype.kind in {"O", "U", "S"}:
            strings = dataset.asstr()[()]
            digest.update("\0".join(map(str, np.asarray(strings).reshape(-1))).encode())
        else:
            digest.update(np.ascontiguousarray(array).tobytes())
    return "sha256:" + digest.hexdigest()


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
        assert handle["/schema/version"].asstr()[()] == "sponge.input.v2"
        identity_uuid = handle["/identity/uuid"].asstr()[()]
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
        assert handle["/schema/version"].asstr()[()] == "sponge.input.v2"
        assert handle["/identity/uuid"].asstr()[()] == identity_uuid
    with h5py.File(paths["restart"], "r") as handle:
        assert handle["/schema/version"].asstr()[()] == "sponge.input.v2"
        assert handle["/identity/uuid"].asstr()[()] == identity_uuid
        assert handle["/particles/all/position/value"].shape == (1, molecule.atom_count, 3)
        assert handle["/particles/all/box/edges/value"].shape == (1, 3, 3)
        assert tuple(handle["/h5md"].attrs["version"]) == (1, 1)
        assert handle["/particles/all/position/value"].attrs["unit"] == "Angstrom"
        assert handle["/particles/all/position/step"].id == handle["/particles/all/step"].id


def test_save_sponge_input_wrapper_selects_raw_or_bundle_format(tmp_path):
    molecule = _peptide()

    assert Xponge.save_sponge_input(molecule, "default", tmp_path) is molecule
    assert (tmp_path / "default_coordinate.txt").is_file()
    assert not (tmp_path / "default_topology.spgt.h5").exists()

    assert Xponge.save_sponge_input(molecule, "raw", tmp_path, format="raw") is molecule
    assert (tmp_path / "raw_coordinate.txt").is_file()

    assert Xponge.save_sponge_input_raw(molecule, "raw_direct", tmp_path) is molecule
    assert (tmp_path / "raw_direct_coordinate.txt").is_file()

    assert Xponge.save_sponge_input(molecule, "bundle", tmp_path, format="bundle") is molecule
    assert (tmp_path / "bundle_topology.spgt.h5").is_file()
    assert (tmp_path / "bundle_protocol.spgp.h5").is_file()
    assert (tmp_path / "bundle_restart.spgr.h5").is_file()


def test_save_sponge_input_bundle_does_not_patch_existing_pdb(tmp_path, monkeypatch):
    molecule = _peptide()
    original_pdb = (
        b"HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        b"CONECT    1    2\n"
        b"END\n"
    )
    pdb_path = tmp_path / "system.pdb"
    pdb_path.write_bytes(original_pdb)
    monkeypatch.chdir(tmp_path)

    Xponge.Save_SPONGE_Input(molecule, "system", tmp_path, format="bundle")

    assert (tmp_path / "system_topology.spgt.h5").is_file()
    assert (tmp_path / "system_protocol.spgp.h5").is_file()
    assert (tmp_path / "system_restart.spgr.h5").is_file()
    assert pdb_path.read_bytes() == original_pdb


def test_save_sponge_input_wrapper_rejects_unknown_format(tmp_path):
    with pytest.raises(ValueError, match="must be 'raw' or 'bundle'"):
        Xponge.save_sponge_input(_peptide(), "system", tmp_path, format="future")


def test_native_bundle_hashes_use_xponge_dataset_sha256_semantics(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()
    Xponge.save_sponge_input_bundle(molecule, "hashed", tmp_path)

    with h5py.File(tmp_path / "hashed_topology.spgt.h5", "r") as topology:
        topology_hash = topology["/topology/topology_hash"].asstr()[()]
        assert topology_hash == _content_hash(topology, "topology.spgt.h5")
        assert topology["/topology/atom_order_hash"].asstr()[()] == _content_hash(
            topology, "topology.spgt.h5", ("/atoms/", "/residues/")
        )
        assert topology["/topology/forcefield_hash"].asstr()[()] == _content_hash(
            topology, "topology.spgt.h5", ("/forcefield/", "/manybody/", "/qc/")
        )
    with h5py.File(tmp_path / "hashed_protocol.spgp.h5", "r") as protocol:
        protocol_hash = protocol["/identity/content_hash"].asstr()[()]
        assert protocol_hash == _content_hash(protocol, "protocol.spgp.h5")
        assert protocol["/protocol/topology_compatibility/topology_hash"].asstr()[()] == topology_hash


def test_native_bundle_topology_hash_changes_for_same_size_content(tmp_path):
    h5py = pytest.importorskip("h5py")
    first = _peptide()
    second = _peptide()
    second.atoms[0].charge += 0.125
    Xponge.save_sponge_input_bundle(first, "first", tmp_path)
    Xponge.save_sponge_input_bundle(second, "second", tmp_path)

    with h5py.File(tmp_path / "first_topology.spgt.h5", "r") as lhs, h5py.File(
        tmp_path / "second_topology.spgt.h5", "r"
    ) as rhs:
        assert lhs["/topology/atom_count"][()] == rhs["/topology/atom_count"][()]
        assert lhs["/topology/topology_hash"].asstr()[()] != rhs["/topology/topology_hash"].asstr()[()]
        assert lhs["/topology/atom_order_hash"].asstr()[()] != rhs["/topology/atom_order_hash"].asstr()[()]
        assert lhs["/topology/forcefield_hash"].asstr()[()] == rhs["/topology/forcefield_hash"].asstr()[()]


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


def test_native_bundle_saver_preserves_supplied_residue_state(tmp_path):
    h5py = pytest.importorskip("h5py")
    parent = _peptide()
    parent.set_box_padding(5.0, center=False)
    residue = parent.residues[0]
    atom = residue.atoms[0]
    atom.x = 12.25
    atom.y = -3.5
    atom.z = 7.75
    atom.charge = -0.4321
    atom.bad_coordinate = True

    molecule = Xponge.save_sponge_input_bundle(residue, "residue", tmp_path)

    assert molecule.residue_count == 1
    assert molecule.atom_count == residue.atom_count
    assert molecule.atoms[0].charge == pytest.approx(-0.4321)
    assert molecule.atoms[0].bad_coordinate is True
    with h5py.File(tmp_path / "residue_topology.spgt.h5", "r") as handle:
        assert handle["/atoms/charge"][0] == pytest.approx(-0.4321 * 18.2223)
    with h5py.File(tmp_path / "residue_restart.spgr.h5", "r") as handle:
        assert handle["/particles/all/position/value"][0, 0].tolist() == pytest.approx(
            [12.25, -3.5, 7.75]
        )


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


def test_native_bundle_saver_materializes_soft_bonds(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()
    molecule.add_bond_soft(1, 0, 12.5, 1.25, 1)

    Xponge.save_sponge_input_bundle(molecule, "soft_bond", tmp_path)

    with h5py.File(tmp_path / "soft_bond_topology.spgt.h5", "r") as handle:
        assert handle["/forcefield/bond_soft/atoms"][:].tolist() == [[0, 1]]
        assert handle["/forcefield/bond_soft/k"][:].tolist() == pytest.approx([12.5])
        assert handle["/forcefield/bond_soft/r0"][:].tolist() == pytest.approx([1.25])
        assert handle["/forcefield/bond_soft/from_a_or_b"][:].tolist() == [1]


def test_native_bundle_saver_materializes_ryckaert_bellemans_listed_force(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _peptide()
    molecule.add_ryckaert_bellemans(4, 3, 2, 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

    Xponge.save_sponge_input_bundle(molecule, "rb", tmp_path)

    root = "/forcefield/custom_force/listed"
    data = f"{root}/data/Ryckaert_Bellemans"
    with h5py.File(tmp_path / "rb_topology.spgt.h5", "r") as handle:
        assert handle[f"{root}/name"].asstr()[:].tolist() == ["Ryckaert_Bellemans"]
        assert handle[f"{root}/parameters/offset"][:].tolist() == [0, 10]
        assert int(handle[f"{root}/count"][()]) == 1
        assert int(handle[f"{data}/item_count"][()]) == 1
        assert handle[f"{data}/parameter/is_int"].dtype == np.dtype(bool)
        assert handle[f"{data}/parameter/int_value"][0, :4].tolist() == [1, 2, 3, 4]
        assert handle[f"{data}/parameter/float_value"][0, 4:].tolist() == pytest.approx(
            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        )
        assert handle["/topology/topology_hash"].asstr()[()] == _content_hash(
            handle, "topology.spgt.h5"
        )


def test_validate_bundle_parity_hash_normalizes_boolean_datasets(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_path = tmp_path / "logical.h5"
    with h5py.File(bundle_path, "w") as handle:
        handle.create_dataset("/logical", data=np.asarray([False, True], dtype=bool))

    helper = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts" / "validate_bundle_parity.py")
    )["_content_hash"]
    with h5py.File(bundle_path, "r") as handle:
        assert helper(handle, "topology.spgt.h5") == _content_hash(
            handle, "topology.spgt.h5"
        )


def test_native_bundle_saver_rejects_minimum_bonded_parameters(tmp_path):
    from XpongeCPP.forcefield.special import min as min_helper

    molecule = _peptide()
    min_helper.save_min_bonded_parameters()
    try:
        with pytest.raises(ValueError, match="minimum-bonded parameters"):
            Xponge.save_sponge_input(
                molecule, "minimum", tmp_path, format="bundle"
            )
    finally:
        min_helper.do_not_save_min_bonded_parameters()

    assert not list(tmp_path.glob("minimum_*"))


def test_native_bundle_saver_rejects_user_defined_listed_forces(tmp_path):
    molecule = _peptide()
    molecule.add_listed_force_definition(
        "[[[ User_Defined ]]]\nfloat k\nE = k;\n[[[ END ]]]\n"
    )

    with pytest.raises(ValueError, match="user-defined listed forces"):
        Xponge.save_sponge_input_bundle(molecule, "custom", tmp_path)

    assert not list(tmp_path.glob("custom_*"))


def test_native_bundle_saver_rejects_stillinger_weber_parameters(tmp_path):
    molecule = _peptide()
    molecule.add_sw_type(
        "S-S", 1.1, 2.2, 3.3, 4.0, 5.0, 6.6, 7.7, 8.8, 0.0, 0.0
    )

    with pytest.raises(ValueError, match="Stillinger-Weber"):
        Xponge.save_sponge_input_bundle(molecule, "sw", tmp_path)

    assert not list(tmp_path.glob("sw_*"))


def test_native_bundle_saver_rejects_edip_parameters(tmp_path):
    molecule = _peptide()
    molecule.add_edip_type(
        "E-E",
        1.1,
        2.2,
        3.3,
        4.4,
        5.5,
        6.6,
        7.7,
        8.8,
        0.0,
        0.0,
        0.0,
        12.12,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )

    with pytest.raises(ValueError, match="EDIP"):
        Xponge.save_sponge_input_bundle(molecule, "edip", tmp_path)

    assert not list(tmp_path.glob("edip_*"))


def test_native_bundle_saver_rejects_escaping_prefix(tmp_path):
    with pytest.raises(ValueError, match="escapes output directory"):
        Xponge.save_sponge_input_bundle(_peptide(), "../escape", tmp_path)


@pytest.mark.parametrize("prefix", [".", "subdir/.."])
def test_native_bundle_saver_rejects_prefix_resolving_to_output_root(
    tmp_path, prefix
):
    output_root = tmp_path / "bundle-root"

    with pytest.raises(ValueError, match="prefix"):
        Xponge.save_sponge_input_bundle(_peptide(), prefix, output_root)

    assert not list(tmp_path.glob("bundle-root_*"))


@pytest.mark.parametrize(
    "angles, message",
    [
        ([90.0, 90.0, 0.0], "singular"),
        ([10.0, 10.0, 170.0], "invalid"),
    ],
)
def test_native_bundle_saver_rejects_invalid_box_geometry(
    tmp_path, angles, message
):
    molecule = _peptide()
    molecule.set_box_padding(5.0)
    molecule.box_angle = angles

    with pytest.raises(ValueError, match=message):
        Xponge.save_sponge_input_bundle(molecule, "invalid_box", tmp_path)

    assert not list(tmp_path.glob("invalid_box_*"))
