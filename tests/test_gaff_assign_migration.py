import io
import json
import re
import site
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
import XpongeCPP as Xponge
from conftest import original_xponge_repo


REPO_ROOT = Path(__file__).resolve().parents[1]
XPONGE_GAFF = original_xponge_repo() / "Xponge" / "forcefield" / "amber" / "gaff.py"
GAFF_100_DIR = REPO_ROOT / "tests" / "data" / "gaff_assign_100"


def _prepare_largest_connected_mol2(source, destination):
    text = Path(source).read_text()
    sections = []
    current_name = None
    current_lines = []
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            if current_name is not None:
                sections.append((current_name, current_lines))
            current_name = line
            current_lines = []
        else:
            current_lines.append(line)
    if current_name is not None:
        sections.append((current_name, current_lines))

    by_name = {name: lines for name, lines in sections}
    atom_lines = [line for line in by_name["@<TRIPOS>ATOM"] if line.strip()]
    bond_lines = [line for line in by_name.get("@<TRIPOS>BOND", []) if line.strip()]
    atom_ids = [int(line.split()[0]) for line in atom_lines]
    adjacency = {atom_id: set() for atom_id in atom_ids}
    parsed_bonds = []
    for line in bond_lines:
        words = line.split()
        atom1, atom2 = int(words[1]), int(words[2])
        parsed_bonds.append((atom1, atom2, words))
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)

    components = []
    remaining = set(atom_ids)
    while remaining:
        root = min(remaining)
        stack = [root]
        component = set()
        while stack:
            atom = stack.pop()
            if atom in component:
                continue
            component.add(atom)
            stack.extend(adjacency[atom] - component)
        remaining -= component
        components.append(component)
    largest = min(components, key=lambda component: (-len(component), min(component)))
    if len(largest) == len(atom_ids):
        return Path(source), False

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    remap = {old: new for new, old in enumerate(sorted(largest), start=1)}
    filtered_atoms = []
    for line in atom_lines:
        words = line.split()
        old_id = int(words[0])
        if old_id not in largest:
            continue
        words[0] = str(remap[old_id])
        filtered_atoms.append(" ".join(words))
    filtered_bonds = []
    for atom1, atom2, words in parsed_bonds:
        if atom1 not in largest or atom2 not in largest:
            continue
        words[0] = str(len(filtered_bonds) + 1)
        words[1] = str(remap[atom1])
        words[2] = str(remap[atom2])
        filtered_bonds.append(" ".join(words))

    molecule_lines = list(by_name["@<TRIPOS>MOLECULE"])
    counts = molecule_lines[1].split()
    counts[0] = str(len(filtered_atoms))
    counts[1] = str(len(filtered_bonds))
    molecule_lines[1] = " ".join(counts)
    replacements = {
        "@<TRIPOS>MOLECULE": molecule_lines,
        "@<TRIPOS>ATOM": filtered_atoms,
        "@<TRIPOS>BOND": filtered_bonds,
    }
    unity_lines = by_name.get("@<TRIPOS>UNITY_ATOM_ATTR")
    if unity_lines is not None:
        filtered_unity = []
        index = 0
        while index < len(unity_lines):
            if not unity_lines[index].strip():
                index += 1
                continue
            words = unity_lines[index].split()
            atom_id, attribute_count = int(words[0]), int(words[1])
            attributes = unity_lines[index + 1:index + 1 + attribute_count]
            if atom_id in largest:
                filtered_unity.append(f"{remap[atom_id]} {attribute_count}")
                filtered_unity.extend(attributes)
            index += 1 + attribute_count
        replacements["@<TRIPOS>UNITY_ATOM_ATTR"] = filtered_unity

    output = []
    for name, lines in sections:
        output.append(name)
        output.extend(replacements.get(name, lines))
    destination.write_text("\n".join(output) + "\n")
    return destination, True


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


def test_gaff_assign_100_matches_current_original_xponge(tmp_path):
    manifest_path = GAFF_100_DIR / "manifest.json"
    if not manifest_path.exists():
        pytest.skip("run benchmarks/generate_gaff_assign_100_baseline.py to create the 100-molecule baseline")

    manifest = json.loads(manifest_path.read_text())
    entries = manifest.get("entries", [])
    assert len(entries) >= 100
    entries = entries[:100]
    prepared = [
        _prepare_largest_connected_mol2(
            GAFF_100_DIR / entry["input_mol2"],
            tmp_path / "prepared" / Path(entry["input_mol2"]).name,
        )
        for entry in entries
    ]
    mol2_paths = [str(path) for path, _ in prepared]
    assert sum(was_prepared for _, was_prepared in prepared) == 8
    reference_path = tmp_path / "xponge-current-gaff-reference.json"
    site_packages = Path(site.getsitepackages()[0])
    script = textwrap.dedent(
        f"""
        import json
        import sys

        sys.path.insert(0, {str(original_xponge_repo())!r})
        sys.path.append({str(site_packages)!r})
        import Xponge
        import Xponge.forcefield.amber.gaff  # noqa: F401

        results = []
        for path in {mol2_paths!r}:
            assignment = Xponge.get_assignment_from_mol2(path, total_charge="sum")
            assignment.determine_atom_type("gaff")
            results.append([
                getattr(assignment.atom_types[index], "name", str(assignment.atom_types[index]))
                for index in range(len(assignment.atoms))
            ])
        with open({str(reference_path)!r}, "w", encoding="utf-8") as handle:
            json.dump(results, handle)
        """
    )
    result = subprocess.run(
        [sys.executable, "-S", "-c", script],
        cwd=original_xponge_repo(),
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        pytest.fail(f"current Xponge GAFF reference failed: {result.stderr[-2000:]}")
    references = json.loads(reference_path.read_text())

    mismatches = []
    for entry, mol2_path, expected in zip(entries, mol2_paths, references):
        assignment = Xponge.get_assignment_from_mol2(mol2_path, total_charge="sum")
        assignment.determine_atom_type("gaff")
        if assignment.atom_types != expected:
            mismatches.append(
                {
                    "source_id": entry.get("source_id", entry.get("cid")),
                    "expected": expected,
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
    assert any(line.startswith("ATOM") for line in pdb_path.read_text().splitlines())

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
