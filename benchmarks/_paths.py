from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA_1KV2_DIR = REPO_ROOT / "tests" / "data" / "1kv2"
PDB_1KV2_H = TEST_DATA_1KV2_DIR / "1KV2_H.pdb"


def original_xponge_repo() -> Path:
    configured = os.environ.get("XPONGE_REFERENCE_REPO")
    if configured:
        return Path(configured)
    candidates = [
        REPO_ROOT.parent / "Xponge-origin",
        REPO_ROOT.parent / "Xponge",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]
