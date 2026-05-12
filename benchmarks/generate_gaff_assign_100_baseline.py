#!/usr/bin/env python3
"""Generate the fixed GAFF Assign parity manifest with original Xponge."""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import re
import signal
import sys
import urllib.request
from pathlib import Path


def _load_paths_module():
    module_path = Path(__file__).resolve().with_name("_paths.py")
    spec = importlib.util.spec_from_file_location("_xpongecpp_bench_paths", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


_paths = _load_paths_module()
REPO_ROOT = _paths.REPO_ROOT
original_xponge_repo = _paths.original_xponge_repo


DEFAULT_XPONGE_REPO = original_xponge_repo()
DEFAULT_CIDS = REPO_ROOT / "tests" / "data" / "gaff_assign_100" / "candidate_cids.txt"
DEFAULT_OUT = REPO_ROOT / "tests" / "data" / "gaff_assign_100"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
CHEMBL_LATEST = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"


def read_cids(path: Path) -> list[int]:
    cids: list[int] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            cids.append(int(line))
    return cids


def import_xponge(xponge_repo: Path, source: str):
    sys.path.insert(0, str(xponge_repo))
    if source == "pubchem":
        try:
            import pubchempy  # noqa: F401  # pylint: disable=unused-import,import-outside-toplevel
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "original Xponge get_assignment_from_pubchem requires pubchempy; "
                "run with `uv run --with pubchempy python benchmarks/generate_gaff_assign_100_baseline.py --source pubchem`"
            ) from exc
    if source == "chembl":
        try:
            import rdkit  # noqa: F401  # pylint: disable=unused-import,import-outside-toplevel
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "original Xponge get_assignment_from_smiles requires rdkit for ChEMBL SMILES; "
                "run with `uv run --with rdkit python benchmarks/generate_gaff_assign_100_baseline.py`"
            ) from exc
    import Xponge  # pylint: disable=import-error,import-outside-toplevel
    import Xponge.forcefield.amber.gaff  # pylint: disable=unused-import,import-error,import-outside-toplevel
    from Xponge.assign import get_assignment_from_pubchem  # pylint: disable=import-error,import-outside-toplevel
    from Xponge.assign import get_assignment_from_smiles  # pylint: disable=import-error,import-outside-toplevel

    return Xponge, get_assignment_from_pubchem, get_assignment_from_smiles


def assignment_from_record(record: dict, get_assignment_from_pubchem, get_assignment_from_smiles):
    if record["source"] == "pubchem":
        return get_assignment_from_pubchem(record["id"], "cid")
    if record["source"] == "chembl":
        return get_assignment_from_smiles(record["smiles"])
    raise ValueError(f"unknown source: {record['source']}")


def molecule_entry(record: dict, get_assignment_from_pubchem, get_assignment_from_smiles,
                   output_dir: Path, input_dir: Path, baseline_mode: str) -> dict:
    source = assignment_from_record(record, get_assignment_from_pubchem, get_assignment_from_smiles)
    source.name = record["label"]
    input_path = input_dir / f"{record['label'].lower()}.mol2"
    source.save_as_mol2(str(input_path), atomtype="sybyl")

    if baseline_mode == "mol2":
        from Xponge.assign import get_assignment_from_mol2  # pylint: disable=import-error,import-outside-toplevel

        typed = get_assignment_from_mol2(str(input_path), total_charge="sum")
    else:
        typed = assignment_from_record(record, get_assignment_from_pubchem, get_assignment_from_smiles)
        typed.name = record["label"]
    typed.determine_atom_type("gaff")

    bonds = []
    for i, atom_bonds in typed.bonds.items():
        for j, order in atom_bonds.items():
            if i < j:
                bonds.append([i, j, order])

    return {
        "source": record["source"],
        "source_id": record["id"],
        "smiles": record.get("smiles"),
        "input_mol2": str(input_path.relative_to(output_dir)),
        "atom_count": typed.atom_numbers,
        "bond_count": len(bonds),
        "baseline_mode": baseline_mode,
        "xponge_gaff_atom_types": [typed.atom_types[i].name for i in range(typed.atom_numbers)],
        "elements": list(typed.atoms),
        "bonds": bonds,
    }


def pubchem_records(candidate_cids: Path) -> list[dict]:
    return [{"source": "pubchem", "id": cid, "label": f"CID{cid}"} for cid in read_cids(candidate_cids)]


def open_url(url: str, timeout_sec: int):
    request = urllib.request.Request(url, headers={"User-Agent": "XpongeCPP GAFF baseline generator"})
    return urllib.request.urlopen(request, timeout=timeout_sec)


def discover_latest_chemreps_url(timeout_sec: int) -> str:
    with open_url(CHEMBL_LATEST, timeout_sec) as handle:
        text = handle.read().decode("utf-8", errors="replace")
    matches = re.findall(r'href="([^"]*chemreps\.txt\.gz)"', text)
    if not matches:
        raise RuntimeError("could not find ChEMBL chemreps txt.gz in latest FTP directory")
    return CHEMBL_LATEST + matches[0]


def chembl_records_from_ftp(limit: int, timeout_sec: int, chemreps_url: str | None) -> list[dict]:
    records: list[dict] = []
    url = chemreps_url or discover_latest_chemreps_url(timeout_sec)
    with open_url(url, timeout_sec) as response:
        with gzip.GzipFile(fileobj=response) as gzip_file:
            header = gzip_file.readline().decode("utf-8", errors="replace").strip().split("\t")
            columns = {name.lower(): index for index, name in enumerate(header)}
            id_column = columns["chembl_id"] if "chembl_id" in columns else columns.get("chemblid")
            smiles_column = (
                columns["canonical_smiles"] if "canonical_smiles" in columns else columns.get("canonicalsmiles")
            )
            if id_column is None or smiles_column is None:
                raise RuntimeError(f"unexpected ChEMBL chemreps header: {header}")
            for raw_line in gzip_file:
                line = raw_line.decode("utf-8", errors="replace").rstrip("\n")
                words = line.split("\t")
                if len(words) <= max(id_column, smiles_column):
                    continue
                chembl_id = words[id_column].strip()
                smiles = words[smiles_column].strip()
                if not chembl_id or not smiles:
                    continue
                records.append({"source": "chembl", "id": chembl_id, "label": chembl_id, "smiles": smiles})
                if len(records) >= limit:
                    return records
    return records


def chembl_records_from_api(limit: int, max_pages: int, page_size: int, timeout_sec: int) -> list[dict]:
    records: list[dict] = []
    seen = set()
    offset = 0
    for _page in range(max_pages):
        query = f"limit={page_size}&offset={offset}&molecule_properties__mw_freebase__lte=250&structure_type=MOL"
        url = CHEMBL_API + "?" + query
        with open_url(url, timeout_sec) as handle:
            payload = json.load(handle)
        for molecule in payload.get("molecules", []):
            chembl_id = molecule.get("molecule_chembl_id")
            structures = molecule.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles")
            if not chembl_id or not smiles or chembl_id in seen:
                continue
            seen.add(chembl_id)
            records.append({"source": "chembl", "id": chembl_id, "label": chembl_id, "smiles": smiles})
            if len(records) >= limit:
                return records
        page_meta = payload.get("page_meta") or {}
        next_url = page_meta.get("next")
        if not next_url:
            break
        offset += page_size
    return records


def run_with_timeout(timeout_sec: int, func, *args):
    if timeout_sec <= 0:
        return func(*args)

    def handler(_signum, _frame):
        raise TimeoutError(f"timed out after {timeout_sec} seconds")

    previous = signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout_sec)
    try:
        return func(*args)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def write_manifest(output_dir: Path, limit: int, entries: list[dict], failures: list[dict]) -> None:
    manifest = {
        "source": "original Xponge Assign + determine_atom_type('gaff')",
        "limit": limit,
        "entries": entries,
        "failures": failures,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xponge-repo", type=Path, default=DEFAULT_XPONGE_REPO)
    parser.add_argument("--source", choices=["chembl", "pubchem"], default="chembl")
    parser.add_argument("--candidate-cids", type=Path, default=DEFAULT_CIDS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--timeout-sec", type=int, default=30)
    parser.add_argument("--chembl-pages", type=int, default=20)
    parser.add_argument("--chembl-page-size", type=int, default=50)
    parser.add_argument("--chembl-mode", choices=["ftp", "api"], default="ftp")
    parser.add_argument("--chembl-chemreps-url", default="")
    parser.add_argument("--candidate-multiplier", type=int, default=10)
    parser.add_argument("--baseline-mode", choices=["direct", "mol2"], default="direct")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    input_dir = output_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    _xponge, get_assignment_from_pubchem, get_assignment_from_smiles = import_xponge(args.xponge_repo.resolve(), args.source)
    entries = []
    failures = []
    if args.source == "chembl":
        candidate_limit = args.limit * args.candidate_multiplier
        if args.chembl_mode == "ftp":
            records = chembl_records_from_ftp(candidate_limit, args.timeout_sec, args.chembl_chemreps_url or None)
        else:
            records = chembl_records_from_api(candidate_limit, args.chembl_pages, args.chembl_page_size, args.timeout_sec)
    else:
        records = pubchem_records(args.candidate_cids)
    for record in records:
        if len(entries) >= args.limit:
            break
        try:
            entries.append(
                run_with_timeout(
                    args.timeout_sec,
                    molecule_entry,
                    record,
                    get_assignment_from_pubchem,
                    get_assignment_from_smiles,
                    output_dir,
                    input_dir,
                    args.baseline_mode,
                )
            )
            print(f"[ok] {record['label']} ({len(entries)}/{args.limit})", flush=True)
        except Exception as exc:  # noqa: BLE001 - baseline generation records original-Xponge failures.
            failures.append({"source": record["source"], "source_id": record["id"], "error": repr(exc)})
            print(f"[skip] {record['label']}: {exc!r}", flush=True)
        write_manifest(output_dir, args.limit, entries, failures)

    if len(entries) < args.limit:
        print(f"only generated {len(entries)} successful entries; need {args.limit}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
