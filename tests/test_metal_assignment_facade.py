from dataclasses import replace
from io import StringIO

import pytest
import XpongeCPP as Xponge
from XpongeCPP.metal_assignment import (
    AtomParameterUpdate,
    BondParameter,
    ChargeUpdate,
    ElectronicState,
    LJParameter,
    MetalAssignmentValidationError,
    MetalSite,
    apply_metal_assignment,
    build_metal_parameter_overlay,
    molecule_input_hash,
    molecule_topology_hash,
    prepare_metal_assignment,
)


METAL_PDB_TEXT = """\
HETATM    1 ZN   MZN A   1       0.000   0.000   0.000  1.00  0.00          Zn
HETATM    2 N1   LIG A   2       2.100   0.000   0.000  1.00  0.00           N
END
"""


def _load_metal_site():
    return Xponge.load_pdb(StringIO(METAL_PDB_TEXT))


def _prepare(molecule):
    return prepare_metal_assignment(
        molecule,
        metal_sites=(MetalSite(atom_id=0, element="Zn", formal_charge=2),),
        electronic_state=ElectronicState(total_charge=0, spin_multiplicity=1),
        coordination_edges=((0, 1),),
        charge_updates=(
            ChargeUpdate(atom_id=0, charge=1.25),
            ChargeUpdate(atom_id=1, charge=-1.25),
        ),
    )


def test_prepare_is_side_effect_free_and_hash_closed():
    molecule = _load_metal_site()
    input_hash = molecule_input_hash(molecule)
    topology_hash = molecule_topology_hash(molecule)
    original_charges = [atom.charge for atom in molecule.atoms]

    plan = _prepare(molecule)

    assert molecule_input_hash(molecule) == input_hash
    assert molecule_topology_hash(molecule) == topology_hash
    assert molecule.residue_links == []
    assert [atom.charge for atom in molecule.atoms] == pytest.approx(
        original_charges
    )
    assert plan.input_hash == input_hash
    assert plan.topology_hash == topology_hash
    assert plan.plan_hash == plan.computed_hash()
    plan.validate_hash()


def test_apply_returns_modified_copy_without_mutating_parent():
    molecule = _load_metal_site()
    parent_hash = molecule_input_hash(molecule)

    result = apply_metal_assignment(molecule, _prepare(molecule))

    assert result.molecule is not molecule
    assert molecule_input_hash(molecule) == parent_hash
    assert molecule.residue_links == []
    assert result.molecule.residue_links == [[0, 1]]
    assert [atom.charge for atom in result.molecule.atoms] == pytest.approx(
        [1.25, -1.25]
    )
    assert result.applied_charge_atom_ids == (0, 1)
    assert result.applied_coordination_edges == ((0, 1),)
    assert result.result_input_hash == molecule_input_hash(result.molecule)
    assert result.result_topology_hash == molecule_topology_hash(result.molecule)
    assert result.inplace is False
    assert result.result_hash == result.computed_hash()
    result.validate_hash()
    assert "returned_copy" in result.application_audit
    assert any(item.startswith("request:") for item in result.provenance)


def test_inplace_apply_publishes_complete_state_once():
    molecule = _load_metal_site()

    result = apply_metal_assignment(molecule, _prepare(molecule), inplace=True)

    assert result.molecule is molecule
    assert molecule.residue_links == [[0, 1]]
    assert [atom.charge for atom in molecule.atoms] == pytest.approx([1.25, -1.25])
    assert result.inplace is True


def test_stale_parent_is_rejected_without_partial_mutation():
    molecule = _load_metal_site()
    plan = _prepare(molecule)
    molecule.atoms[1].charge = 0.5
    state_before_apply = molecule_input_hash(molecule)

    with pytest.raises(MetalAssignmentValidationError) as exc:
        apply_metal_assignment(molecule, plan, inplace=True)

    assert exc.value.code == "stale_parent_molecule"
    assert molecule_input_hash(molecule) == state_before_apply
    assert molecule.residue_links == []


def test_tampered_plan_is_rejected_before_mutation():
    molecule = _load_metal_site()
    plan = _prepare(molecule)
    tampered = replace(plan, link_overlay=())
    parent_hash = molecule_input_hash(molecule)

    with pytest.raises(MetalAssignmentValidationError) as exc:
        apply_metal_assignment(molecule, tampered, inplace=True)

    assert exc.value.code == "stale_plan_hash"
    assert molecule_input_hash(molecule) == parent_hash
    assert molecule.residue_links == []


@pytest.mark.parametrize(
    ("kwargs", "code"),
    [
        ({"metal_sites": ()}, "missing_metal_sites"),
        ({"coordination_edges": ()}, "missing_coordination_edges"),
        (
            {"metal_sites": (MetalSite(0, "Fe", 2),)},
            "metal_element_mismatch",
        ),
        (
            {"electronic_state": ElectronicState(0, 0)},
            "invalid_spin_multiplicity",
        ),
    ],
)
def test_prepare_rejects_incomplete_or_inconsistent_inputs(kwargs, code):
    molecule = _load_metal_site()
    arguments = {
        "metal_sites": (MetalSite(0, "Zn", 2),),
        "electronic_state": ElectronicState(0, 1),
        "coordination_edges": ((0, 1),),
    }
    arguments.update(kwargs)

    with pytest.raises(MetalAssignmentValidationError) as exc:
        prepare_metal_assignment(molecule, **arguments)

    assert exc.value.code == code


def test_legacy_namespace_exports_molecule_first_facade():
    from Xponge.metal_assignment import prepare_metal_assignment as legacy_prepare

    assert legacy_prepare is prepare_metal_assignment


def test_parameter_overlay_is_hash_closed_and_applied_locally():
    molecule = _load_metal_site()
    original_type = molecule.atoms[0].type
    original_mass = molecule.atoms[0].mass
    overlay = build_metal_parameter_overlay(
        molecule,
        atom_parameters=(
            AtomParameterUpdate(
                atom_id=0,
                atom_type="metal_Zn_test",
                mass=65.4,
                source="fixture",
            ),
        ),
        bond_parameters=(
            BondParameter(
                atom_ids=(0, 1),
                force_constant=123.5,
                equilibrium_length=2.1,
                source="fixture",
            ),
        ),
        parameter_source="fixture",
    )
    plan = prepare_metal_assignment(
        molecule,
        metal_sites=(MetalSite(0, "Zn", 2),),
        electronic_state=ElectronicState(0, 1),
        coordination_edges=((0, 1),),
        parameter_overlay=overlay,
    )

    result = apply_metal_assignment(molecule, plan)

    assert overlay.overlay_hash == overlay.computed_hash()
    assert result.molecule.atoms[0].type == "metal_Zn_test"
    assert result.molecule.atoms[0].mass == pytest.approx(65.4)
    assert molecule.atoms[0].type == original_type
    assert molecule.atoms[0].mass == pytest.approx(original_mass)


def test_parameter_overlay_reaches_raw_and_bundle_bond_savers(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _load_metal_site()
    overlay = build_metal_parameter_overlay(
        molecule,
        atom_parameters=(
            AtomParameterUpdate(0, "metal_Zn_export_test", 65.4, "fixture"),
            AtomParameterUpdate(1, "donor_N_export_test", 14.01, "fixture"),
        ),
        bond_parameters=(
            BondParameter((0, 1), 123.5, 2.1, "fixture"),
        ),
        lj_parameters=(
            LJParameter(
                "metal_Zn_export_test",
                "metal_Zn_export_test",
                0.0125,
                1.1,
                "fixture",
            ),
            LJParameter(
                "donor_N_export_test",
                "donor_N_export_test",
                0.025,
                1.5,
                "fixture",
            ),
        ),
        parameter_source="fixture",
    )
    plan = prepare_metal_assignment(
        molecule,
        metal_sites=(MetalSite(0, "Zn", 2),),
        electronic_state=ElectronicState(0, 1),
        coordination_edges=((0, 1),),
        parameter_overlay=overlay,
    )
    applied = apply_metal_assignment(molecule, plan).molecule

    Xponge.Save_SPONGE_Input(applied, "metal", dirname=str(tmp_path / "raw"))
    Xponge.save_sponge_input_bundle(applied, "metal", tmp_path / "bundle")

    raw_row = (tmp_path / "raw" / "metal_bond.txt").read_text().splitlines()[1]
    assert raw_row == "0 1 123.500000 2.100000"
    assert len(applied.lj_parameter_overrides) == 2
    with h5py.File(
        tmp_path / "bundle" / "metal_topology.spgt.h5", "r"
    ) as handle:
        assert handle["/forcefield/bond/atoms"][:].tolist() == [[0, 1]]
        assert handle["/forcefield/bond/k"][:].tolist() == pytest.approx([123.5])
        assert handle["/forcefield/bond/r0"][:].tolist() == pytest.approx([2.1])


def test_tampered_parameter_overlay_is_rejected_before_mutation():
    molecule = _load_metal_site()
    overlay = build_metal_parameter_overlay(
        molecule,
        bond_parameters=(BondParameter((0, 1), 100.0, 2.0, "fixture"),),
    )
    tampered = replace(
        overlay,
        bond_parameters=(BondParameter((0, 1), 999.0, 2.0, "fixture"),),
    )

    with pytest.raises(MetalAssignmentValidationError) as exc:
        prepare_metal_assignment(
            molecule,
            metal_sites=(MetalSite(0, "Zn", 2),),
            electronic_state=ElectronicState(0, 1),
            coordination_edges=((0, 1),),
            parameter_overlay=tampered,
        )

    assert exc.value.code == "stale_parameter_overlay_hash"
    assert molecule.residue_links == []
