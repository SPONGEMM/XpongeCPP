"""MDAnalysis topology support for CIF/mmCIF paired with H5MD."""

from __future__ import annotations

import hashlib
import json
import os
import warnings
from typing import Any

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis.core.topology import Topology
from MDAnalysis.lib.util import openany
from MDAnalysis.topology.base import TopologyReaderBase

from ..io_bundle.errors import BundleTopologyError, BundleTrajectoryError


_MMCIF_MISSING_VALUES = {"", ".", "?"}


def _mmcif_clean(value):
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in _MMCIF_MISSING_VALUES else text


def _mmcif_first(row, names):
    for name in names:
        value = _mmcif_clean(row.get(name.lower()))
        if value:
            return value
    return ""


def _mmcif_parse_int(value):
    value = _mmcif_clean(value)
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _mmcif_parse_float(value):
    value = _mmcif_clean(value)
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _mmcif_chain_id(value):
    value = _mmcif_clean(value)
    return value[0] if value else " "


def _mmcif_insertion_code(value):
    value = _mmcif_clean(value)
    return value[0] if value else " "


def _mmcif_tokenize(text):
    tokens = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            block = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                block.append(lines[i])
                i += 1
            if i >= len(lines):
                raise ValueError("unterminated mmCIF text field")
            tokens.append("\n".join(block))
            i += 1
            continue
        j = 0
        current = []
        quote = None
        while j < len(line):
            char = line[j]
            if quote:
                if char == quote and (j + 1 == len(line) or line[j + 1].isspace()):
                    tokens.append("".join(current))
                    current = []
                    quote = None
                else:
                    current.append(char)
                j += 1
                continue
            if char == "#":
                break
            if char.isspace():
                if current:
                    tokens.append("".join(current))
                    current = []
                j += 1
                continue
            if char in ("'", '"'):
                if current:
                    current.append(char)
                else:
                    quote = char
                j += 1
                continue
            current.append(char)
            j += 1
        if quote:
            raise ValueError("unterminated mmCIF quoted value")
        if current:
            tokens.append("".join(current))
        i += 1
    return tokens


def _mmcif_parse(file):
    text = file.read()
    tokens = _mmcif_tokenize(text)
    data = {}
    i = 0
    while i < len(tokens):
        token = tokens[i]
        lower = token.lower()
        if lower.startswith("data_") or lower.startswith("save_"):
            i += 1
            continue
        if lower == "loop_":
            i += 1
            tags = []
            while i < len(tokens) and tokens[i].startswith("_"):
                tags.append(tokens[i].lower())
                i += 1
            if not tags:
                raise ValueError("mmCIF loop without tags")
            values = []
            while i < len(tokens):
                next_lower = tokens[i].lower()
                at_row_boundary = len(values) % len(tags) == 0
                if at_row_boundary and (
                    next_lower == "loop_" or next_lower.startswith("data_") or
                    next_lower.startswith("save_") or tokens[i].startswith("_")
                ):
                    break
                values.append(tokens[i])
                i += 1
            if len(values) % len(tags) != 0:
                raise ValueError("mmCIF loop row has incomplete values")
            for tag_index, tag in enumerate(tags):
                data[tag] = values[tag_index::len(tags)]
            continue
        if token.startswith("_"):
            if i + 1 >= len(tokens):
                raise ValueError(f"mmCIF tag without value: {token}")
            data[token.lower()] = [tokens[i + 1]]
            i += 2
            continue
        i += 1
    return data


def _mmcif_rows(data, category):
    prefix = "_" + category.lower() + "."
    tags = [tag for tag in data if tag.startswith(prefix)]
    if not tags:
        return []
    count = max(len(data[tag]) for tag in tags)
    rows = []
    for index in range(count):
        row = {}
        for tag in tags:
            values = data[tag]
            row[tag] = values[index] if index < len(values) else ""
        rows.append(row)
    return rows



def _mokda_mapping_digest(mapping_rows: list[dict[str, Any]]) -> str:
    """Return Mokda's canonical digest for ordered simulation atom mappings."""

    fields = (
        "simulation_index",
        "external_id",
        "canonical_atom_id",
        "simulation_residue_index",
        "simulation_residue_id",
    )
    payload = {
        "atom_mapping": [
            {field: row[field] for field in fields} for row in mapping_rows
        ]
    }
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


class CIFTopologyParser(TopologyReaderBase):
    """Read mmCIF ``_atom_site`` metadata into an MDAnalysis topology.

    The parser follows Xponge's mmCIF conventions: author atom/residue/chain
    identifiers take precedence, label identifiers fill missing author fields,
    and alternate location ``A`` is selected by default.  Explicit PDB segment
    IDs are retained as MDAnalysis ``segids``; otherwise ``label_asym_id`` is
    used as the segment ID while the author chain remains available as
    ``chainID``.
    """

    format = "XPONGE_CIF"

    @staticmethod
    def _format_hint(thing: Any) -> bool:
        if not isinstance(thing, (str, os.PathLike)):
            return False
        filename = os.fspath(thing)
        lower_filename = filename.lower()
        if lower_filename.endswith(
            (".cif", ".mmcif", ".mcif", ".cif.gz", ".mmcif.gz", ".mcif.gz")
        ):
            return True
        try:
            with openany(filename, "rt") as handle:
                return "_atom_site." in handle.read(65536).lower()
        except (OSError, TypeError, ValueError):
            return False

    def parse(self, **kwargs) -> Topology:
        with openany(self.filename, "rt") as handle:
            data = _mmcif_parse(handle)
        atom_rows = _mmcif_rows(data, "atom_site")
        if not atom_rows:
            raise ValueError("mmCIF file does not contain _atom_site records")

        def first(row: dict[str, str], *tags: str) -> str:
            return _mmcif_first(row, list(tags))

        def parse_formal_charge(value: str) -> int | None:
            cleaned = _mmcif_clean(value)
            charge = _mmcif_parse_int(cleaned)
            if charge is None and len(cleaned) > 1 and cleaned[-1] in "+-":
                magnitude = _mmcif_parse_int(cleaned[:-1])
                if magnitude is not None:
                    charge = magnitude if cleaned[-1] == "+" else -magnitude
            return charge

        model_values = {
            first(row, "_atom_site.pdbx_pdb_model_num") or "1" for row in atom_rows
        }
        model_id = kwargs.get("model_id")
        if model_id is None and len(model_values) > 1:
            raise ValueError("mmCIF contains multiple models; pass model_id explicitly")
        selected_model = str(model_id) if model_id is not None else next(iter(model_values))
        if selected_model not in model_values:
            raise ValueError(f"mmCIF does not contain model {selected_model!r}")

        altloc = kwargs.get("altloc", "A")
        selected_rows = []
        for row in atom_rows:
            row_model = first(row, "_atom_site.pdbx_pdb_model_num") or "1"
            if row_model != selected_model:
                continue
            row_altloc = first(row, "_atom_site.label_alt_id")
            if altloc is not None and row_altloc and row_altloc != str(altloc):
                continue
            selected_rows.append(row)
        if not selected_rows:
            raise ValueError("mmCIF model/alternate-location selection produced no atoms")

        atom_ids: list[int] = []
        atom_names: list[str] = []
        atom_types: list[str] = []
        elements: list[str] = []
        residue_names: list[str] = []
        residue_ids: list[int] = []
        residue_numbers: list[int] = []
        residue_index_by_atom: list[int] = []
        residue_segment_index: list[int] = []
        segment_ids: list[str] = []
        chain_ids: list[str] = []
        insertion_codes: list[str] = []
        alternate_locations: list[str] = []
        occupancies: list[float | None] = []
        temperature_factors: list[float | None] = []
        formal_charges: list[int | None] = []
        record_types: list[str] = []
        residue_components: list[str] = []
        residue_insertion_codes: list[str] = []
        residue_atom_names: list[dict[str, int]] = []
        atom_by_site_id: dict[str, int] = {}
        atom_by_label_key: dict[tuple[str, ...], int] = {}
        atom_by_auth_key: dict[tuple[str, ...], int] = {}
        segment_index_by_id: dict[str, int] = {}
        residue_index_by_key: dict[tuple[str, ...], int] = {}
        atom_id_set: set[int] = set()

        def identity_key(asym_id, seq_id, comp_id, atom_id, insertion_code):
            return (
                _mmcif_clean(asym_id),
                _mmcif_clean(seq_id),
                _mmcif_clean(comp_id).upper(),
                _mmcif_clean(atom_id).upper(),
                _mmcif_clean(insertion_code),
            )

        for row in selected_rows:
            auth_asym = first(row, "_atom_site.auth_asym_id")
            label_asym = first(row, "_atom_site.label_asym_id")
            chain_id = auth_asym or label_asym
            segid = (
                first(
                    row,
                    "_atom_site.pdbx_pdb_segment_id",
                    "_atom_site.seg_id",
                    "_atom_site.segid",
                )
                or label_asym
                or auth_asym
                or "SYSTEM"
            )
            if segid not in segment_index_by_id:
                segment_index_by_id[segid] = len(segment_ids)
                segment_ids.append(segid)
            segment_index = segment_index_by_id[segid]

            auth_seq = first(row, "_atom_site.auth_seq_id")
            label_seq = first(row, "_atom_site.label_seq_id")
            auth_resid = _mmcif_parse_int(auth_seq)
            label_resid = _mmcif_parse_int(label_seq)
            if auth_resid is None:
                auth_resid = label_resid
            if label_resid is None:
                label_resid = auth_resid
            auth_comp = first(row, "_atom_site.auth_comp_id")
            label_comp = first(row, "_atom_site.label_comp_id")
            residue_name = auth_comp or label_comp or "SYSTEM"
            insertion_code = first(row, "_atom_site.pdbx_pdb_ins_code")
            residue_key = (
                segid,
                chain_id,
                label_asym,
                auth_seq,
                label_seq,
                insertion_code,
                residue_name,
            )
            if residue_key not in residue_index_by_key:
                residue_index = len(residue_names)
                residue_index_by_key[residue_key] = residue_index
                resid = auth_resid if auth_resid is not None else residue_index + 1
                resnum = label_resid if label_resid is not None else resid
                residue_names.append(residue_name)
                residue_ids.append(resid)
                residue_numbers.append(resnum)
                residue_insertion_codes.append(insertion_code)
                residue_segment_index.append(segment_index)
                residue_components.append(label_comp or auth_comp or residue_name)
                residue_atom_names.append({})
            residue_index = residue_index_by_key[residue_key]
            residue_index_by_atom.append(residue_index)

            atom_index = len(atom_names)
            raw_atom_id = first(row, "_atom_site.id")
            parsed_atom_id = _mmcif_parse_int(raw_atom_id)
            if parsed_atom_id is None or parsed_atom_id in atom_id_set:
                parsed_atom_id = atom_index + 1
                while parsed_atom_id in atom_id_set:
                    parsed_atom_id += 1
            atom_id_set.add(parsed_atom_id)
            atom_ids.append(parsed_atom_id)

            auth_atom = first(row, "_atom_site.auth_atom_id")
            label_atom = first(row, "_atom_site.label_atom_id")
            atom_name = auth_atom or label_atom or f"A{atom_index + 1}"
            atom_names.append(atom_name)
            type_symbol = first(row, "_atom_site.type_symbol")
            element = type_symbol[:1].upper() + type_symbol[1:].lower() if type_symbol else ""
            elements.append(element)
            atom_types.append(element or atom_name)
            chain_ids.append(chain_id)
            insertion_codes.append(insertion_code)
            alternate_locations.append(first(row, "_atom_site.label_alt_id"))
            occupancies.append(
                _mmcif_parse_float(first(row, "_atom_site.occupancy"))
            )
            temperature_factors.append(
                _mmcif_parse_float(first(row, "_atom_site.b_iso_or_equiv"))
            )
            formal_charge = parse_formal_charge(
                first(row, "_atom_site.pdbx_formal_charge")
            )
            formal_charges.append(formal_charge)
            record_types.append(first(row, "_atom_site.group_pdb").upper() or "ATOM")

            label_seq = first(row, "_atom_site.label_seq_id")
            label_comp = first(row, "_atom_site.label_comp_id")
            label_key = identity_key(label_asym, label_seq, label_comp, label_atom, insertion_code)
            auth_key = identity_key(
                auth_asym or label_asym,
                auth_seq or label_seq,
                auth_comp or label_comp,
                auth_atom or label_atom,
                insertion_code,
            )
            if any(label_key):
                atom_by_label_key[label_key] = atom_index
            if any(auth_key):
                atom_by_auth_key[auth_key] = atom_index
            if raw_atom_id:
                atom_by_site_id[raw_atom_id] = atom_index
            for alias in (auth_atom, label_atom):
                if alias:
                    residue_atom_names[residue_index][alias.upper()] = atom_index

        atom_count = len(atom_names)
        residue_count = len(residue_names)
        segment_count = len(segment_ids)

        trajectory_mapping = []
        mapping_rows = _mmcif_rows(data, "mokda_trajectory_mapping")
        if mapping_rows:
            mapping_by_index: dict[int, dict[str, Any]] = {}
            mapping_hashes: set[str] = set()
            canonical_atom_ids: set[int] = set()
            external_ids: set[str] = set()
            simulation_residue_id_by_index: dict[int, str] = {}
            simulation_residue_index_by_id: dict[str, int] = {}
            for row in mapping_rows:
                simulation_index = _mmcif_parse_int(
                    first(row, "_mokda_trajectory_mapping.simulation_index")
                )
                atom_site_id = first(row, "_mokda_trajectory_mapping.atom_site_id")
                canonical_atom_id = _mmcif_parse_int(
                    first(row, "_mokda_trajectory_mapping.canonical_atom_id")
                )
                external_id = first(row, "_mokda_trajectory_mapping.external_id")
                simulation_residue_index = _mmcif_parse_int(
                    first(row, "_mokda_trajectory_mapping.simulation_residue_index")
                )
                simulation_residue_id = first(
                    row, "_mokda_trajectory_mapping.simulation_residue_id"
                )
                mapping_hash = first(row, "_mokda_trajectory_mapping.mapping_hash")
                if (
                    simulation_index is None
                    or canonical_atom_id is None
                    or simulation_residue_index is None
                    or atom_site_id not in atom_by_site_id
                    or not external_id
                    or not simulation_residue_id
                    or canonical_atom_id <= 0
                    or simulation_residue_index < 0
                ):
                    raise BundleTopologyError(
                        "invalid _mokda_trajectory_mapping row"
                    )
                if simulation_index in mapping_by_index:
                    raise BundleTopologyError(
                        "duplicate _mokda_trajectory_mapping.simulation_index"
                    )
                if canonical_atom_id in canonical_atom_ids:
                    raise BundleTopologyError(
                        "duplicate _mokda_trajectory_mapping.canonical_atom_id"
                    )
                if external_id in external_ids:
                    raise BundleTopologyError(
                        "duplicate _mokda_trajectory_mapping.external_id"
                    )
                canonical_atom_ids.add(canonical_atom_id)
                external_ids.add(external_id)
                previous_residue_id = simulation_residue_id_by_index.setdefault(
                    simulation_residue_index, simulation_residue_id
                )
                previous_residue_index = simulation_residue_index_by_id.setdefault(
                    simulation_residue_id, simulation_residue_index
                )
                if (
                    previous_residue_id != simulation_residue_id
                    or previous_residue_index != simulation_residue_index
                ):
                    raise BundleTopologyError(
                        "conflicting _mokda_trajectory_mapping simulation residue identity"
                    )
                mapping_by_index[simulation_index] = {
                    "simulation_index": simulation_index,
                    "atom_site_id": atom_site_id,
                    "canonical_atom_id": canonical_atom_id,
                    "external_id": external_id,
                    "simulation_residue_index": simulation_residue_index,
                    "simulation_residue_id": simulation_residue_id,
                    "mapping_hash": mapping_hash,
                }
                mapping_hashes.add(mapping_hash)

            expected_indices = list(range(len(selected_rows)))
            if sorted(mapping_by_index) != expected_indices:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping simulation indices must cover "
                    "the selected CIF atoms from 0 to n_atoms - 1"
                )
            trajectory_mapping = [
                mapping_by_index[index] for index in expected_indices
            ]
            mapped_site_ids = [row["atom_site_id"] for row in trajectory_mapping]
            selected_site_ids = [first(row, "_atom_site.id") for row in selected_rows]
            if len(set(selected_site_ids)) != len(selected_site_ids):
                raise BundleTopologyError(
                    "_atom_site.id values must be unique when a Mokda trajectory mapping is present"
                )
            if mapped_site_ids != selected_site_ids:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping does not match _atom_site row order"
                )
            if len(mapping_hashes) != 1:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping rows have inconsistent mapping hashes"
                )
            mapping_hash = next(iter(mapping_hashes))
            if (
                len(mapping_hash) != 64
                or mapping_hash != mapping_hash.lower()
                or any(character not in "0123456789abcdef" for character in mapping_hash)
            ):
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping.mapping_hash must be lowercase SHA-256"
                )
            if _mokda_mapping_digest(trajectory_mapping) != mapping_hash:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping rows do not match mapping_hash"
                )

        partial_charges: list[float | None] = [None] * atom_count
        for row in _mmcif_rows(data, "mokda_charge_state"):
            site_id = first(row, "_mokda_charge_state.atom_site_id")
            atom_index = atom_by_site_id.get(site_id)
            if atom_index is None:
                # A multi-model CIF may contain charge rows for atoms outside
                # the model selected for this topology.
                continue
            partial_charges[atom_index] = _mmcif_parse_float(
                first(row, "_mokda_charge_state.partial_charge")
            )
            if formal_charges[atom_index] is None:
                formal_charges[atom_index] = parse_formal_charge(
                    first(row, "_mokda_charge_state.formal_charge")
                )

        attrs: list[Any] = [
            topologyattrs.Atomids(np.asarray(atom_ids, dtype=np.int64)),
            topologyattrs.Atomnames(np.asarray(atom_names, dtype=object)),
            topologyattrs.Atomtypes(np.asarray(atom_types, dtype=object)),
            topologyattrs.Resids(np.asarray(residue_ids, dtype=np.int64)),
            topologyattrs.Resnums(np.asarray(residue_numbers, dtype=np.int64)),
            topologyattrs.Resnames(np.asarray(residue_names, dtype=object)),
            topologyattrs.Segids(np.asarray(segment_ids, dtype=object)),
            topologyattrs.ChainIDs(np.asarray(chain_ids, dtype=object)),
            topologyattrs.ICodes(np.asarray(residue_insertion_codes, dtype=object)),
            topologyattrs.RecordTypes(np.asarray(record_types, dtype=object)),
        ]
        if any(elements):
            attrs.append(topologyattrs.Elements(np.asarray(elements, dtype=object)))
        if any(alternate_locations):
            attrs.append(topologyattrs.AltLocs(np.asarray(alternate_locations, dtype=object)))
        if any(value is not None for value in occupancies):
            attrs.append(
                topologyattrs.Occupancies(
                    np.asarray([np.nan if value is None else value for value in occupancies], dtype=np.float64)
                )
            )
        if any(value is not None for value in temperature_factors):
            attrs.append(
                topologyattrs.Tempfactors(
                    np.asarray(
                        [np.nan if value is None else value for value in temperature_factors],
                        dtype=np.float64,
                    )
                )
            )
        if any(value is not None for value in formal_charges):
            attrs.append(
                topologyattrs.FormalCharges(
                    np.asarray([0 if value is None else value for value in formal_charges], dtype=np.int64)
                )
            )
        if any(value is not None for value in partial_charges):
            attrs.append(
                topologyattrs.Charges(
                    np.asarray(
                        [np.nan if value is None else value for value in partial_charges],
                        dtype=np.float64,
                    )
                )
            )

        bonds: set[tuple[int, int]] = set()
        bond_orders: dict[tuple[int, int], float] = {}
        bond_types: dict[tuple[int, int], str] = {}

        def parse_bond_order(value: str) -> float | None:
            token = _mmcif_clean(value).upper()
            if token in {"SING", "SINGLE"}:
                return 1.0
            if token in {"DOUB", "DOUBLE"}:
                return 2.0
            if token in {"TRIP", "TRIPLE"}:
                return 3.0
            if token in {"AROM", "AROMATIC", "DELO"}:
                return 1.5
            return _mmcif_parse_float(token)

        def add_bond(
            atom1: int | None,
            atom2: int | None,
            *,
            order: float | None = None,
            bond_type: str = "",
        ) -> None:
            if atom1 is None or atom2 is None or atom1 == atom2:
                return
            bond = tuple(sorted((int(atom1), int(atom2))))
            bonds.add(bond)
            if order is not None:
                bond_orders[bond] = order
            if bond_type:
                bond_types[bond] = bond_type

        def remove_bond(atom1: int | None, atom2: int | None) -> None:
            if atom1 is None or atom2 is None or atom1 == atom2:
                return
            bond = tuple(sorted((int(atom1), int(atom2))))
            bonds.discard(bond)
            bond_orders.pop(bond, None)
            bond_types.pop(bond, None)

        for row in _mmcif_rows(data, "chem_comp_bond"):
            comp_id = first(row, "_chem_comp_bond.comp_id").upper()
            atom1_name = first(row, "_chem_comp_bond.atom_id_1").upper()
            atom2_name = first(row, "_chem_comp_bond.atom_id_2").upper()
            if not comp_id or not atom1_name or not atom2_name:
                continue
            for residue_index, residue_comp in enumerate(residue_components):
                if residue_comp.upper() == comp_id:
                    add_bond(
                        residue_atom_names[residue_index].get(atom1_name),
                        residue_atom_names[residue_index].get(atom2_name),
                        order=parse_bond_order(
                            first(row, "_chem_comp_bond.value_order")
                        ),
                    )

        def resolve_struct_conn_atom(row: dict[str, str], partner: str) -> int | None:
            label_key = identity_key(
                first(row, f"_struct_conn.{partner}_label_asym_id"),
                first(row, f"_struct_conn.{partner}_label_seq_id"),
                first(row, f"_struct_conn.{partner}_label_comp_id"),
                first(row, f"_struct_conn.{partner}_label_atom_id"),
                first(row, f"_struct_conn.pdbx_{partner}_pdb_ins_code"),
            )
            auth_key = identity_key(
                first(row, f"_struct_conn.{partner}_auth_asym_id", f"_struct_conn.{partner}_label_asym_id"),
                first(row, f"_struct_conn.{partner}_auth_seq_id", f"_struct_conn.{partner}_label_seq_id"),
                first(row, f"_struct_conn.{partner}_auth_comp_id", f"_struct_conn.{partner}_label_comp_id"),
                first(row, f"_struct_conn.{partner}_auth_atom_id", f"_struct_conn.{partner}_label_atom_id"),
                first(row, f"_struct_conn.pdbx_{partner}_pdb_ins_code"),
            )
            return atom_by_auth_key.get(auth_key, atom_by_label_key.get(label_key))

        for row in _mmcif_rows(data, "struct_conn"):
            add_bond(
                resolve_struct_conn_atom(row, "ptnr1"),
                resolve_struct_conn_atom(row, "ptnr2"),
                bond_type=first(row, "_struct_conn.conn_type_id"),
            )
        for row in _mmcif_rows(data, "mokda_bond_semantic"):
            add_bond(
                atom_by_site_id.get(first(row, "_mokda_bond_semantic.atom_site_id_1")),
                atom_by_site_id.get(first(row, "_mokda_bond_semantic.atom_site_id_2")),
                order=_mmcif_parse_float(
                    first(row, "_mokda_bond_semantic.order")
                ),
                bond_type=first(row, "_mokda_bond_semantic.bond_type"),
            )
        for row in _mmcif_rows(data, "mokda_edit_operation"):
            operation = first(row, "_mokda_edit_operation.operation_type").lower()
            atom1 = atom_by_site_id.get(first(row, "_mokda_edit_operation.atom_site_id_1"))
            atom2 = atom_by_site_id.get(first(row, "_mokda_edit_operation.atom_site_id_2"))
            if atom1 is None or atom2 is None or atom1 == atom2:
                continue
            bond = tuple(sorted((int(atom1), int(atom2))))
            if operation == "delete_bond":
                remove_bond(atom1, atom2)
            elif operation in {"add_bond", "create_bond", "update_bond"}:
                add_bond(
                    atom1,
                    atom2,
                    order=_mmcif_parse_float(
                        first(row, "_mokda_edit_operation.bond_order")
                    ),
                    bond_type=first(row, "_mokda_edit_operation.bond_type"),
                )
        if bonds:
            ordered_bonds = sorted(bonds)
            bond_kwargs = {}
            if bond_orders:
                bond_kwargs["order"] = np.asarray(
                    [bond_orders.get(bond, np.nan) for bond in ordered_bonds],
                    dtype=np.float64,
                )
            if bond_types:
                bond_kwargs["types"] = np.asarray(
                    [bond_types.get(bond, "") for bond in ordered_bonds],
                    dtype=object,
                )
            attrs.append(topologyattrs.Bonds(ordered_bonds, **bond_kwargs))

        topology = Topology(
            atom_count,
            residue_count,
            segment_count,
            attrs,
            np.asarray(residue_index_by_atom, dtype=np.int32),
            np.asarray(residue_segment_index, dtype=np.int32),
        )
        topology.cif_metadata = {
            # Keep every parsed CIF column, including Mokda-specific categories
            # and auth/label identifiers that have no one-to-one MDA attribute.
            "raw": data,
            "selected_model": selected_model,
            "selected_atom_site_ids": tuple(
                first(row, "_atom_site.id") for row in selected_rows
            ),
            "mokda_trajectory_mapping": trajectory_mapping,
        }
        return topology


def _mokda_mapping_file_for_cif(topology) -> str | None:
    """Find Mokda's standard mapping sidecar beside an ordered topology CIF."""

    if not isinstance(topology, (str, os.PathLike)):
        return None
    filename = os.fspath(topology)
    suffix = "_trajectory_topology.cif"
    basename = os.path.basename(filename)
    if not basename.endswith(suffix):
        return None
    prefix = basename[: -len(suffix)]
    return os.path.join(
        os.path.dirname(filename), f"{prefix}_atom_order_mapping.json"
    )


def _verify_mokda_mapping_file(mapping_file, cif_mapping):
    """Compare Mokda's JSON mapping sidecar with its CIF mapping category."""

    if not cif_mapping:
        raise BundleTopologyError(
            "a Mokda atom-order mapping file is present, but the CIF has no "
            "_mokda_trajectory_mapping category"
        )
    try:
        with open(mapping_file, "rt", encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BundleTopologyError(
            f"could not read Mokda atom-order mapping file {mapping_file!r}: {exc}"
        ) from exc
    if (
        not isinstance(document, dict)
        or document.get("schema") != "sponge-atom-order-mapping"
    ):
        raise BundleTopologyError(
            f"{mapping_file!r} is not a sponge-atom-order-mapping document"
        )
    mapping_atoms = document.get("atoms")
    if not isinstance(mapping_atoms, list):
        raise BundleTopologyError(
            f"{mapping_file!r} does not contain an atoms list"
        )

    fields = (
        "simulation_index",
        "external_id",
        "canonical_atom_id",
        "simulation_residue_index",
        "simulation_residue_id",
    )
    try:
        sidecar_mapping = [
            {field: row[field] for field in fields}
            for row in mapping_atoms
            if isinstance(row, dict)
        ]
    except KeyError as exc:
        raise BundleTopologyError(
            f"{mapping_file!r} has an incomplete atom mapping row"
        ) from exc
    if len(sidecar_mapping) != len(mapping_atoms):
        raise BundleTopologyError(
            f"{mapping_file!r} contains a non-object atom mapping row"
        )
    cif_identity_mapping = [
        {field: row[field] for field in fields} for row in cif_mapping
    ]
    cif_mapping_hash = cif_mapping[0]["mapping_hash"]
    if document.get("mapping_hash") != cif_mapping_hash:
        raise BundleTopologyError(
            "Mokda mapping file hash does not match the CIF trajectory mapping"
        )
    if sidecar_mapping != cif_identity_mapping:
        raise BundleTopologyError(
            "Mokda mapping file atoms do not match the CIF trajectory mapping"
        )
    return {
        "path": os.fspath(mapping_file),
        "verified": True,
        "mapping_hash": cif_mapping_hash,
        "atom_count": len(cif_identity_mapping),
    }


def load_cif_h5md_universe(
    topology,
    trajectory,
    *,
    mapping_file=None,
    model_id: str | int | None = None,
    altloc: str | None = "A",
    particle_stream: str | None = None,
    walker: int | None = None,
    strict: bool = True,
    allow_incomplete: bool = False,
    convert_units: bool = True,
    **universe_kwargs,
):
    """Load a CIF topology and H5MD trajectory into an MDAnalysis Universe.

    Mokda's companion atom-order JSON is checked when it exists. A missing
    sidecar is allowed; the CIF mapping category, when present, is still
    validated against atom-site row order and its embedded digest.
    """
    from .bundle_mdanalysis import SpongeH5MDReader

    if not CIFTopologyParser._format_hint(topology):
        raise ValueError(f"{topology!r} is not a recognized CIF/mmCIF topology path")
    try:
        trajectory_atom_count = SpongeH5MDReader.parse_n_atoms(
            trajectory,
            particle_stream=particle_stream,
            walker=walker,
        )
    except (OSError, KeyError, ValueError) as exc:
        raise BundleTrajectoryError(
            f"{trajectory!r} is not a recognized H5MD trajectory"
        ) from exc

    cif_topology = CIFTopologyParser(topology).parse(model_id=model_id, altloc=altloc)
    if mapping_file is None:
        mapping_file = _mokda_mapping_file_for_cif(topology)
    if mapping_file is not None and os.path.isfile(mapping_file):
        cif_topology.cif_metadata["mapping_file_validation"] = (
            _verify_mokda_mapping_file(
                mapping_file,
                cif_topology.cif_metadata["mokda_trajectory_mapping"],
            )
        )
    else:
        cif_topology.cif_metadata["mapping_file_validation"] = {
            "available": False,
            "verified": False,
        }

    if trajectory_atom_count != cif_topology.n_atoms:
        raise BundleTrajectoryError(
            f"{trajectory!r} has {trajectory_atom_count} atoms, CIF topology has "
            f"{cif_topology.n_atoms}"
        )
    if cif_topology.cif_metadata["mokda_trajectory_mapping"]:
        mapping_file_verified = cif_topology.cif_metadata["mapping_file_validation"][
            "verified"
        ]
        mapping_status = (
            "the CIF simulation mapping matches its companion mapping file"
            if mapping_file_verified
            else "the CIF simulation mapping matches its atom-site row order and digest"
        )
        warning = (
            f"CIF/H5MD compatibility: {mapping_status}, and the H5MD atom count "
            "matches; this H5MD file has no per-atom IDs for an independent "
            "coordinate-order check."
        )
    else:
        warning = (
            "CIF/H5MD compatibility is checked by atom count only; ensure the CIF "
            "atom order matches the trajectory atom order."
        )
    warnings.warn(warning, RuntimeWarning, stacklevel=2)

    universe = mda.Universe(
        cif_topology,
        trajectory,
        format=SpongeH5MDReader,
        particle_stream=particle_stream,
        walker=walker,
        strict=strict,
        allow_incomplete=allow_incomplete,
        convert_units=convert_units,
        **universe_kwargs,
    )
    universe.cif_metadata = cif_topology.cif_metadata
    return universe


__all__ = [
    "CIFTopologyParser",
    "load_cif_h5md_universe",
]
