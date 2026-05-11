import json
import re
from pathlib import Path

import pytest
import XpongeCPP as Xponge


REPO_ROOT = Path(__file__).resolve().parents[1]
XPONGE_GAFF = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge/Xponge/forcefield/amber/gaff.py")
GAFF_100_DIR = REPO_ROOT / "tests" / "data" / "gaff_assign_100"


def _original_gaff_rule_names():
    if not XPONGE_GAFF.exists():
        pytest.skip("local Xponge GAFF reference is not available")
    text = XPONGE_GAFF.read_text()
    return re.findall(r'@gaff\.Add_Rule\("([^"]+)"\)', text)


def test_current_gaff_assign_types_are_original_xponge_rule_names():
    original = set(_original_gaff_rule_names())
    implemented = set(Xponge.implemented_gaff_assign_types())
    assert implemented <= original


def test_gaff_assign_rule_coverage_matches_original_xponge():
    assert Xponge.implemented_gaff_assign_types() == _original_gaff_rule_names()


def test_gaff_assign_100_manifest_matches_original_xponge_baseline():
    manifest_path = GAFF_100_DIR / "manifest.json"
    if not manifest_path.exists():
        pytest.skip("run benchmarks/generate_gaff_assign_100_baseline.py to create the 100-molecule baseline")

    manifest = json.loads(manifest_path.read_text())
    entries = manifest.get("entries", [])
    assert len(entries) >= 100

    mismatches = []
    for entry in entries[:100]:
        mol2_path = GAFF_100_DIR / entry["input_mol2"]
        assignment = Xponge.get_assignment_from_mol2(str(mol2_path), total_charge="sum")
        assignment.determine_atom_type("gaff")
        if assignment.atom_types != entry["xponge_gaff_atom_types"]:
            mismatches.append(
                {
                    "source_id": entry.get("source_id", entry.get("cid")),
                    "expected": entry["xponge_gaff_atom_types"],
                    "actual": assignment.atom_types,
                }
            )
    assert mismatches == []
