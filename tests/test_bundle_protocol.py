from pathlib import Path

import h5py
import pytest

import XpongeCPP as Xponge
from XpongeCPP.io_bundle.errors import BundlePathError, BundleValidationError


def _peptide():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.get_template_molecule("ALA")


def test_protocol_remapping_preserves_input_atom_semantics():
    from dataclasses import replace
    from XpongeCPP.io_bundle.protocol import remap_protocol_atom_order

    original = _full_protocol(4)
    original = replace(original, collective_variables=(replace(
        original.collective_variables[0], reference_coordinates=((1., 2., 3.), (4., 5., 6.)),
    ),))
    mapped = remap_protocol_atom_order(original, (2, 3, 0, 1))
    assert mapped.collective_variables[0].atom_indices == (2, 3)
    assert mapped.collective_variables[0].reference_coordinates == original.collective_variables[0].reference_coordinates
    assert mapped.distance_constraints[0].atoms == ((2, 3),)
    assert mapped.positional_restraints[0].atom_indices == (2,)
    assert mapped.positional_restraints[0].reference_coordinates == tuple(
        original.positional_restraints[0].reference_coordinates[i] for i in (2, 3, 0, 1)
    )
    assert mapped.positional_restraints[0].weight == original.positional_restraints[0].weight
    assert mapped.sits.atom_indices == (2, 3)
    assert original.sits.atom_indices == (0, 1)
    assert mapped.metadynamics == original.metadynamics
    assert mapped.cv_restraints == original.cv_restraints
    assert mapped.steering == original.steering
    assert remap_protocol_atom_order(original, (0, 1, 2, 3)) is original
    prefix = replace(original, sits=replace(original.sits, atom_indices=(), atom_numbers_policy=2))
    assert remap_protocol_atom_order(prefix, (2, 3, 0, 1)).sits.atom_indices == (2, 3)
    assert remap_protocol_atom_order(prefix, (2, 3, 0, 1)).sits.atom_numbers_policy is None
    with pytest.raises(BundleValidationError, match='permutation'):
        remap_protocol_atom_order(original, (0, 0, 2, 3))


def _full_protocol(atom_count):
    reference = tuple((float(index), 0.0, 0.0) for index in range(atom_count))
    return Xponge.SpongeProtocol(
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="distance_cv",
                type="distance",
                atom_indices=(0, 1),
                period=(6.283185,),
                sigma=(0.1,),
            ),
        ),
        distance_constraints=(
            Xponge.ProtocolDistanceConstraints(
                atoms=((0, 1),),
                r0=(1.45,),
            ),
        ),
        positional_restraints=(
            Xponge.ProtocolPositionalRestraint(
                name="position",
                atom_indices=(0,),
                reference_coordinates=reference,
                weight=((1.0, 1.0, 1.0),),
            ),
        ),
        cv_restraints=(
            Xponge.ProtocolCVRestraint(
                name="cv_bias",
                cv_refs=("distance_cv",),
                weight=(2.0,),
                reference=(1.5,),
            ),
        ),
        metadynamics=(
            Xponge.ProtocolMetadynamics(
                name="meta",
                cv_refs=("distance_cv",),
                grid_min=(0.0,),
                grid_max=(4.0,),
                grid_count=(41,),
                hill_height_default=0.5,
            ),
        ),
        steering=Xponge.ProtocolSteering(
            cv_refs=("distance_cv",),
            weight=(0.25,),
        ),
        sits=Xponge.ProtocolSITS(
            mode="observation",
            atom_indices=(0, 1),
        ),
        hard_wall=Xponge.ProtocolHardWall(
            bounds_low=(0.0, None, None),
            bounds_high=(10.0, None, None),
        ),
        soft_walls=(
            Xponge.ProtocolSoftWall(
                name="soft",
                potential="0.5 * x * x",
            ),
        ),
    )


def _text(dataset):
    value = dataset[()]
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def test_bundle_saver_remaps_protocol_with_reordered_atoms(tmp_path):
    from io import StringIO
    from test_xpongecpp_api import CUSTOM_MOL2_TEXT, MOL2_TEXT

    Xponge.register_tip3p()
    molecule = Xponge.load_mol2(StringIO(CUSTOM_MOL2_TEXT)) | Xponge.load_mol2(StringIO(MOL2_TEXT))
    molecule.add_residue_link(molecule.residues[0].name2atom("O1"), molecule.residues[2].name2atom("O"))
    reference = tuple((float(i), 1., 2.) for i in range(molecule.atom_count))
    protocol = Xponge.SpongeProtocol(
        collective_variables=(Xponge.ProtocolCollectiveVariable(
            name="distance", type="distance", atom_indices=(0, 3),
        ),),
        positional_restraints=(Xponge.ProtocolPositionalRestraint(
            name="position", atom_indices=(3,), reference_coordinates=reference, single_weight_default=1.,
        ),),
    )
    _, mapping = Xponge.save_sponge_input_bundle(
        molecule, "reordered", tmp_path, protocol=protocol,
        source_atom_ids=tuple(range(molecule.atom_count)), return_mapping=True,
    )
    order = [int(row["source_atom_id"]) for row in mapping]
    with h5py.File(tmp_path / "reordered_protocol.spgp.h5") as handle:
        assert handle["/cv/distance/atom_indices"][:].tolist() == [order.index(0), order.index(3)]
    with h5py.File(tmp_path / "reordered_restart.spgr.h5") as handle:
        assert handle["/parameters/restart/references/restraint/position/coordinate"][:].tolist() == [list(reference[i]) for i in order]


def test_native_bundle_protocol_writes_all_origin_contract_groups(tmp_path):
    molecule = _peptide()
    protocol = _full_protocol(molecule.atom_count)

    result = Xponge.save_sponge_input_bundle(
        molecule, "protocol", tmp_path, protocol=protocol
    )

    assert result.atom_count == molecule.atom_count
    protocol_path = tmp_path / "protocol_protocol.spgp.h5"
    restart_path = tmp_path / "protocol_restart.spgr.h5"
    with h5py.File(protocol_path, "r") as handle:
        for path in (
            "/cv/distance_cv/type",
            "/constraint/default/pairs/atoms",
            "/restraint/position/type",
            "/restraint/cv_bias/cv_refs",
            "/meta/meta/grid/count",
            "/steer/cv_refs",
            "/sits/method/mode",
            "/wall/hard/bounds_low",
            "/wall/soft/potential",
        ):
            assert path in handle
        assert int(handle["/protocol/cv_count"][()]) == 1
        assert int(handle["/protocol/restraint_count"][()]) == 2
        protocol_hash = _text(handle["/identity/content_hash"])
    with h5py.File(restart_path, "r") as handle:
        assert (
            "/parameters/restart/references/restraint/position/coordinate"
            in handle
        )
        assert _text(handle["/run/producer_protocol_hash"]) == protocol_hash


def test_invalid_protocol_does_not_publish_partial_bundle(tmp_path):
    molecule = _peptide()
    protocol = Xponge.SpongeProtocol(
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="bad",
                type="distance",
                atom_indices=(0, molecule.atom_count),
            ),
        )
    )

    with pytest.raises(BundleValidationError, match="atom indices"):
        Xponge.save_sponge_input_bundle(
            molecule, "invalid", tmp_path, protocol=protocol
        )

    assert not list(tmp_path.glob("invalid_*"))


def test_protocol_writes_typed_virtual_atom_arrays_and_cv_references(tmp_path):
    molecule = _peptide()
    protocol = Xponge.SpongeProtocol(
        virtual_atoms=(
            Xponge.ProtocolVirtualAtom(
                name="center",
                type="center",
                atom_indices=(0, 1),
                weight=(0.25, 0.75),
            ),
        ),
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="distance_cv",
                type="distance",
                atom_refs=("center", 2),
            ),
        ),
    )

    Xponge.save_sponge_input_bundle(molecule, "virtual", tmp_path, protocol=protocol)

    with h5py.File(tmp_path / "virtual_protocol.spgp.h5", "r") as handle:
        assert _text(handle["/cv/virtual_atom/center/type"]) == "center"
        assert handle["/cv/virtual_atom/center/atom_indices"][...].tolist() == [0, 1]
        assert handle["/cv/virtual_atom/center/weight"][...].tolist() == pytest.approx([0.25, 0.75])
        assert handle["/cv/distance_cv/atom_refs"].asstr()[...].tolist() == ["center", "2"]
        assert int(handle["/protocol/cv_count"][()]) == 1


@pytest.mark.parametrize("prefix", ["../escape", ".", "subdir/.."])
def test_protocol_saver_rejects_unsafe_prefixes(tmp_path, prefix):
    with pytest.raises(BundlePathError, match="prefix|escapes"):
        Xponge.save_sponge_input_bundle(_peptide(), prefix, tmp_path)

    assert not any(
        path.name.startswith("escape_") for path in Path(tmp_path).parent.iterdir()
    )
