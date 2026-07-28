from dataclasses import replace
from io import StringIO

import pytest
import XpongeCPP as Xponge
from XpongeCPP.metal_assignment import (
    ChargeUpdate,
    ElectronicState,
    MetalAssignmentValidationError,
    MetalSite,
    apply_metal_assignment,
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

    plan = _prepare(molecule)

    assert molecule_input_hash(molecule) == input_hash
    assert molecule_topology_hash(molecule) == topology_hash
    assert molecule.residue_links == []
    assert [atom.charge for atom in molecule.atoms] == pytest.approx([0.0, 0.0])
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
