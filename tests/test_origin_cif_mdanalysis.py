from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest


mda = pytest.importorskip("MDAnalysis")

import XpongeCPP.analysis as analysis  # noqa: E402
import XpongeCPP.analysis.md_analysis as xmda  # noqa: E402
from XpongeCPP.io_bundle.errors import (  # noqa: E402
    BundleTopologyError,
    BundleTrajectoryError,
)
from XpongeCPP.tools.traj_analysis import load_Sponge_trajectory  # noqa: E402


def _mapping_document():
    atoms = [
        {
            "simulation_index": 0,
            "external_id": "atom:101",
            "canonical_atom_id": 101,
            "simulation_residue_index": 0,
            "simulation_residue_id": "residue:1",
        },
        {
            "simulation_index": 1,
            "external_id": "atom:202",
            "canonical_atom_id": 202,
            "simulation_residue_index": 0,
            "simulation_residue_id": "residue:1",
        },
        {
            "simulation_index": 2,
            "external_id": "atom:303",
            "canonical_atom_id": 303,
            "simulation_residue_index": 1,
            "simulation_residue_id": "residue:2",
        },
    ]
    canonical = json.dumps(
        {"atom_mapping": atoms},
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return {
        "schema": "sponge-atom-order-mapping",
        "schema_version": 1,
        "mapping_hash": hashlib.sha256(canonical).hexdigest(),
        "atoms": atoms,
    }


def _write_cif(path, *, mapping_site_ids=("1", "2", "3"), include_mapping=True):
    mapping = _mapping_document()
    atom_rows = [
        ("ATOM", "1", "C", "C1", "A", "CMP", "LABEL_A", "7", "1", "A",
         "1.234", "2.345", "3.456", "1.0", "20.0", "0", "XP_SEG",
         "CA", "CMP", "CHAIN_A", "100", "1"),
        ("HETATM", "2", "O", "O1", "A", "CMP", "LABEL_A", "7", "1", "A",
         "4.567", "5.678", "6.789", "0.8", "21.0", "0", "XP_SEG",
         "OX", "CMP", "CHAIN_A", "100", "1"),
        ("HETATM", "3", "N", "N1", "A", "CMP2", "LABEL_B", "8", "2", "?",
         "7.890", "8.901", "9.012", "1.0", "22.0", "?", "?",
         "N1", "CMP2", "CHAIN_B", "7", "1"),
        # A second alternate location is retained in raw CIF metadata but
        # excluded from the default altloc A topology.
        ("HETATM", "4", "O", "O1", "B", "CMP", "LABEL_A", "7", "1", "A",
         "4.567", "5.678", "6.789", "0.2", "25.0", "0", "XP_SEG",
         "OX", "CMP", "CHAIN_A", "100", "1"),
    ]
    atom_tags = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "pdbx_PDB_segment_id",
        "auth_atom_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_model_num",
    ]
    sections = [
        "data_trajectory_topology",
        "_entry.id trajectory_topology",
        "#",
        "loop_",
        *[f"_atom_site.{tag}" for tag in atom_tags],
        *[" ".join(row) for row in atom_rows],
        "#",
        "_mokda_atom_modeling.modeling_role ligand",
        "#",
        "loop_",
        "_mokda_charge_state.atom_site_id",
        "_mokda_charge_state.partial_charge",
        "_mokda_charge_state.formal_charge",
        "1 0.25 0",
        "2 -0.25 0",
        "3 -0.75 2+",
        "#",
        "loop_",
        "_chem_comp_bond.comp_id",
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "CMP CA OX DOUB",
        "#",
        "loop_",
        "_struct_conn.conn_type_id",
        "_struct_conn.ptnr1_label_asym_id",
        "_struct_conn.ptnr1_label_seq_id",
        "_struct_conn.ptnr1_label_comp_id",
        "_struct_conn.ptnr1_label_atom_id",
        "_struct_conn.ptnr1_auth_asym_id",
        "_struct_conn.ptnr1_auth_seq_id",
        "_struct_conn.ptnr1_auth_comp_id",
        "_struct_conn.ptnr1_auth_atom_id",
        "_struct_conn.pdbx_ptnr1_pdb_ins_code",
        "_struct_conn.ptnr2_label_asym_id",
        "_struct_conn.ptnr2_label_seq_id",
        "_struct_conn.ptnr2_label_comp_id",
        "_struct_conn.ptnr2_label_atom_id",
        "_struct_conn.ptnr2_auth_asym_id",
        "_struct_conn.ptnr2_auth_seq_id",
        "_struct_conn.ptnr2_auth_comp_id",
        "_struct_conn.ptnr2_auth_atom_id",
        "_struct_conn.pdbx_ptnr2_pdb_ins_code",
        "covale LABEL_A 1 CMP C1 CHAIN_A 100 CMP CA A "
        "LABEL_B 2 CMP2 N1 CHAIN_B 7 CMP2 N1 ?",
        "#",
        "loop_",
        "_mokda_bond_semantic.atom_site_id_1",
        "_mokda_bond_semantic.atom_site_id_2",
        "_mokda_bond_semantic.order",
        "_mokda_bond_semantic.bond_type",
        "2 3 1.0 covalent",
        "#",
        "loop_",
        "_mokda_edit_operation.atom_site_id_1",
        "_mokda_edit_operation.atom_site_id_2",
        "_mokda_edit_operation.operation_type",
        "_mokda_edit_operation.bond_order",
        "_mokda_edit_operation.bond_type",
        "2 3 update_bond 1.5 aromatic",
        "#",
    ]
    if include_mapping:
        sections.extend(
            [
                "loop_",
                "_mokda_trajectory_mapping.simulation_index",
                "_mokda_trajectory_mapping.atom_site_id",
                "_mokda_trajectory_mapping.canonical_atom_id",
                "_mokda_trajectory_mapping.external_id",
                "_mokda_trajectory_mapping.simulation_residue_index",
                "_mokda_trajectory_mapping.simulation_residue_id",
                "_mokda_trajectory_mapping.mapping_hash",
            ]
        )
        for atom, atom_site_id in zip(mapping["atoms"], mapping_site_ids):
            sections.append(
                f"{atom['simulation_index']} {atom_site_id} "
                f"{atom['canonical_atom_id']} {atom['external_id']} "
                f"{atom['simulation_residue_index']} "
                f"{atom['simulation_residue_id']} {mapping['mapping_hash']}"
            )
        sections.append("#")
    path.write_text("\n".join(sections) + "\n", encoding="utf-8")


def _write_h5md(path, *, atom_count=3):
    source = mda.Universe.empty(atom_count, trajectory=True)
    with mda.coordinates.H5MD.H5MDWriter(
        path,
        n_atoms=atom_count,
        positions=True,
        velocities=False,
        forces=False,
        lengthunit="Angstrom",
        timeunit="ps",
    ) as writer:
        for frame in range(2):
            source.atoms.positions = np.asarray(
                [[index + frame, index + 1, index + 2] for index in range(atom_count)],
                dtype=np.float32,
            )
            writer.write(source)


def test_cif_h5md_route_preserves_atom_fields_mapping_and_bonds(tmp_path):
    topology = tmp_path / "system_trajectory_topology.cif"
    trajectory = tmp_path / "trajectory.h5md"
    _write_cif(topology)
    _write_h5md(trajectory)

    with pytest.warns(RuntimeWarning, match="CIF/H5MD compatibility"):
        universe = load_Sponge_trajectory(topology, trajectory, box=None)

    assert universe.trajectory.__class__ is xmda.SpongeH5MDReader
    assert universe.atoms.ids.tolist() == [1, 2, 3]
    assert universe.atoms.names.tolist() == ["CA", "OX", "N1"]
    assert universe.atoms.types.tolist() == ["C", "O", "N"]
    assert universe.atoms.elements.tolist() == ["C", "O", "N"]
    assert universe.atoms.record_types.tolist() == ["ATOM", "HETATM", "HETATM"]
    assert universe.atoms.chainIDs.tolist() == ["CHAIN_A", "CHAIN_A", "CHAIN_B"]
    assert universe.atoms.segids.tolist() == ["XP_SEG", "XP_SEG", "LABEL_B"]
    assert universe.atoms.resnames.tolist() == ["CMP", "CMP", "CMP2"]
    assert universe.atoms.resids.tolist() == [100, 100, 7]
    assert universe.atoms.resnums.tolist() == [1, 1, 2]
    assert universe.atoms.icodes.tolist() == ["A", "A", ""]
    assert universe.atoms.altLocs.tolist() == ["A", "A", "A"]
    assert universe.atoms.occupancies == pytest.approx([1.0, 0.8, 1.0])
    assert universe.atoms.tempfactors == pytest.approx([20.0, 21.0, 22.0])
    assert universe.atoms.formalcharges.tolist() == [0, 0, 2]
    assert universe.atoms.charges == pytest.approx([0.25, -0.25, -0.75])

    assert len(universe.bonds) == 3
    assert [bond.order for bond in universe.bonds] == pytest.approx([2.0, np.nan, 1.5], nan_ok=True)
    assert universe.bonds[0].type == ""
    assert universe.bonds[1].type == "covale"
    assert universe.bonds[2].type == "aromatic"

    assert universe.cif_metadata["raw"]["_atom_site.label_entity_id"] == [
        "7", "7", "8", "7"
    ]
    assert universe.cif_metadata["raw"]["_atom_site.cartn_x"] == [
        "1.234", "4.567", "7.890", "4.567"
    ]
    assert universe.cif_metadata["raw"]["_mokda_atom_modeling.modeling_role"] == [
        "ligand"
    ]
    assert [row["simulation_index"] for row in universe.cif_metadata["mokda_trajectory_mapping"]] == [
        0, 1, 2
    ]
    assert universe.cif_metadata["mapping_file_validation"] == {
        "available": False,
        "verified": False,
    }

    universe.trajectory[1]
    assert universe.atoms.positions[:, 0] == pytest.approx([1.0, 2.0, 3.0])


def test_cif_without_mokda_mapping_or_sidecar_is_accepted(tmp_path):
    topology = tmp_path / "ordinary.cif"
    trajectory = tmp_path / "ordinary.h5md"
    _write_cif(topology, include_mapping=False)
    _write_h5md(trajectory)

    with pytest.warns(RuntimeWarning, match="atom count only"):
        universe = xmda.load_cif_h5md_universe(topology, trajectory)

    assert universe.atoms.n_atoms == 3
    assert universe.cif_metadata["mokda_trajectory_mapping"] == []
    assert universe.cif_metadata["mapping_file_validation"] == {
        "available": False,
        "verified": False,
    }


def test_cif_mokda_mapping_must_match_atom_site_order(tmp_path):
    topology = tmp_path / "system_trajectory_topology.cif"
    _write_cif(topology, mapping_site_ids=("2", "1", "3"))

    with pytest.raises(BundleTopologyError, match="does not match _atom_site row order"):
        xmda.CIFTopologyParser(topology).parse()


def test_cif_h5md_rejects_atom_count_mismatch(tmp_path):
    topology = tmp_path / "system_trajectory_topology.cif"
    trajectory = tmp_path / "trajectory.h5md"
    _write_cif(topology)
    _write_h5md(trajectory, atom_count=4)

    with pytest.raises(BundleTrajectoryError, match="has 4 atoms, CIF topology has 3"):
        xmda.load_cif_h5md_universe(topology, trajectory)


def test_cif_companion_mapping_is_optional_but_checked_when_present(tmp_path):
    topology = tmp_path / "system_trajectory_topology.cif"
    trajectory = tmp_path / "trajectory.h5md"
    mapping_path = tmp_path / "system_atom_order_mapping.json"
    mapping = _mapping_document()
    _write_cif(topology)
    _write_h5md(trajectory)
    mapping_path.write_text(json.dumps(mapping), encoding="utf-8")

    with pytest.warns(RuntimeWarning, match="companion mapping file"):
        universe = xmda.load_cif_h5md_universe(topology, trajectory)

    assert universe.cif_metadata["mapping_file_validation"]["verified"] is True
    assert universe.cif_metadata["mapping_file_validation"]["atom_count"] == 3

    mapping["atoms"][1]["external_id"] = "atom:wrong"
    mapping_path.write_text(json.dumps(mapping), encoding="utf-8")
    with pytest.raises(BundleTopologyError, match="atoms do not match"):
        xmda.load_cif_h5md_universe(topology, trajectory)


def test_cif_reader_is_public_from_analysis_package():
    assert analysis.CIFTopologyParser is xmda.CIFTopologyParser
    assert analysis.load_cif_h5md_universe is xmda.load_cif_h5md_universe


def test_mda_universe_pairs_cif_topology_reader_with_h5md_reader(tmp_path):
    topology = tmp_path / "system_trajectory_topology.cif"
    trajectory = tmp_path / "trajectory.h5md"
    _write_cif(topology)
    _write_h5md(trajectory)

    universe = mda.Universe(
        topology,
        trajectory,
        topology_format=xmda.CIFTopologyParser,
        format=xmda.SpongeH5MDReader,
    )

    assert universe.atoms.ids.tolist() == [1, 2, 3]
    assert universe.atoms.names.tolist() == ["CA", "OX", "N1"]
    assert universe.trajectory.n_atoms == 3
