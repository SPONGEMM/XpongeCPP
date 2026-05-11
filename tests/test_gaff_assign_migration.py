import io
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


def test_assignment_from_mol2_preserves_sybyl_atom_type_details():
    mol2 = """@<TRIPOS>MOLECULE
SYBYL_TYPES
 4 3 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1 C1      0.0000   0.0000   0.0000 C.ar       1 MOL 0.000000
     2 N1      1.3000   0.0000   0.0000 N.pl3      1 MOL 0.000000
     3 O1      2.6000   0.0000   0.0000 O.co2      1 MOL 0.000000
     4 CL1     3.9000   0.0000   0.0000 Cl         1 MOL 0.000000
@<TRIPOS>BOND
     1 1 2 ar
     2 2 3 1
     3 3 4 1
@<TRIPOS>SUBSTRUCTURE
     1 MOL 1
"""

    assignment = Xponge.get_assignment_from_mol2(io.StringIO(mol2), total_charge="sum")

    assert assignment.atoms == ["C", "N", "O", "Cl"]
    assert assignment.element_details == [".ar", ".pl3", ".co2", ""]

    assignment.determine_atom_type("sybyl")
    assert assignment.atom_types == ["C.ar", "N.pl3", "O.co2", "Cl"]


def test_assignment_compatibility_entrypoints_and_writers(tmp_path):
    xyz = """3
water
O 0.0000 0.0000 0.0000
H 0.9572 0.0000 0.0000
H -0.2399 0.9266 0.0000
"""
    assignment = Xponge.get_assignment_from_xyz(io.StringIO(xyz))
    assert assignment.atoms == ["O", "H", "H"]
    assert assignment.atom_count == 3

    mol2_path = tmp_path / "assigned.mol2"
    pdb_path = tmp_path / "assigned.pdb"
    assignment.save_as_mol2(str(mol2_path), residue_name="WAT")
    assignment.save_as_pdb(str(pdb_path), residue_name="WAT")
    assert "@<TRIPOS>ATOM" in mol2_path.read_text()
    assert pdb_path.read_text().startswith("ATOM")

    pdb_assignment = Xponge.get_assignment_from_pdb(io.StringIO(pdb_path.read_text()))
    assert pdb_assignment.atoms == ["O", "H", "H"]

    residue_type = assignment.to_residuetype("WAT")
    residue_assignment = Xponge.get_assignment_from_residuetype(residue_type)
    assert residue_assignment.atoms == ["O", "H", "H"]


def test_smiles_and_pubchem_entrypoints_report_missing_optional_dependencies_clearly(monkeypatch):
    import builtins

    original_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("rdkit") or name.startswith("pubchempy"):
            raise ImportError(name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ImportError, match="RDKit"):
        Xponge.get_assignment_from_smiles("CCO")
    with pytest.raises(ImportError, match="PubChemPy"):
        Xponge.get_assignment_from_pubchem("ethanol")
