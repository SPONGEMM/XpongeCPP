import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA_DIR = REPO_ROOT / "tests" / "data"
DATA_1KV2_DIR = TEST_DATA_DIR / "1kv2"


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


def optional_1kv2_baseline_dir() -> Path | None:
    configured = os.environ.get("XPONGECPP_1KV2_BASELINE_DIR")
    if not configured:
        candidates = [
            TEST_DATA_DIR / "1kv2_xponge_baseline",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None
    return Path(configured)
