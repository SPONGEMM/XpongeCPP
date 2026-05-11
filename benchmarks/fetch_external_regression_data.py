#!/usr/bin/env python3
"""Fetch reproducible external regression inputs for XpongeCPP.

Network access is intentionally kept out of the normal test suite. This script
creates cached fixtures and manifests that tests can consume offline.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import NamedTuple


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RCSB_CACHE = REPO_ROOT / "tests" / "data" / "rcsb_cache"
DEFAULT_CHEMBL_1000 = REPO_ROOT / "tests" / "data" / "gaff_assign_1000"
RCSB_IDS = ("1CRN", "1UBQ", "1AKE", "4HHB", "1BNA")
USER_AGENT = "XpongeCPP external regression data fetcher"


class RcsbEntry(NamedTuple):
    pdb_id: str
    url: str
    cache_name: str


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rcsb-cache-dir", type=Path, default=DEFAULT_RCSB_CACHE)
    parser.add_argument("--chembl-output-dir", type=Path, default=DEFAULT_CHEMBL_1000)
    parser.add_argument("--chembl-limit", type=int, default=1000)
    parser.add_argument("--timeout-sec", type=int, default=30)
    parser.add_argument("--xponge-repo", type=Path, default=REPO_ROOT.parent / "Xponge")
    parser.add_argument("--skip-rcsb", action="store_true")
    parser.add_argument("--skip-chembl", action="store_true")
    parser.add_argument("--baseline-mode", choices=["direct", "mol2"], default="mol2")
    return parser


def planned_rcsb_entries() -> list[RcsbEntry]:
    return [
        RcsbEntry(pdb_id=pdb_id, url=f"https://files.rcsb.org/download/{pdb_id}.pdb", cache_name=f"{pdb_id}.pdb")
        for pdb_id in RCSB_IDS
    ]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cached_file_manifest(entry: RcsbEntry, path: Path, fetched_at: str) -> dict:
    return {
        "pdb_id": entry.pdb_id,
        "url": entry.url,
        "cache_path": path.name,
        "fetched_at": fetched_at,
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def fetch_url(url: str, path: Path, timeout_sec: int) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout_sec) as response:
        path.write_bytes(response.read())


def fetch_rcsb(cache_dir: Path, timeout_sec: int) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    fetched_at = dt.date.today().isoformat()
    manifest = []
    for entry in planned_rcsb_entries():
        path = cache_dir / entry.cache_name
        if not path.exists():
            fetch_url(entry.url, path, timeout_sec)
        manifest.append(cached_file_manifest(entry, path, fetched_at))
    (cache_dir / "manifest.json").write_text(json.dumps({"entries": manifest}, indent=2) + "\n")


def generate_chembl_baseline(output_dir: Path, limit: int, timeout_sec: int, xponge_repo: Path, baseline_mode: str) -> None:
    generator = REPO_ROOT / "benchmarks" / "generate_gaff_assign_100_baseline.py"
    command = [
        sys.executable,
        str(generator),
        "--source",
        "chembl",
        "--limit",
        str(limit),
        "--output-dir",
        str(output_dir),
        "--timeout-sec",
        str(timeout_sec),
        "--xponge-repo",
        str(xponge_repo),
        "--baseline-mode",
        baseline_mode,
    ]
    subprocess.run(command, check=True)


def main() -> int:
    args = build_parser().parse_args()
    if not args.skip_rcsb:
        fetch_rcsb(args.rcsb_cache_dir, args.timeout_sec)
    if not args.skip_chembl:
        generate_chembl_baseline(
            args.chembl_output_dir,
            args.chembl_limit,
            args.timeout_sec,
            args.xponge_repo,
            args.baseline_mode,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
