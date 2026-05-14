import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_8RYK_DIR = REPO_ROOT / "tests" / "data" / "8ryk"
LEGACY_8RYK_ROOT = "/media/yuh/BCDC9249DC91FDB8/Data/Mokda-FEP/8RYK"


def _snapshot_files(path: Path) -> dict[str, str]:
    return {
        item.name: item.read_text()
        for item in sorted(path.iterdir())
        if item.is_file()
    }


def _bond_adjacency_from_text(text: str) -> dict[int, set[int]]:
    adjacency = {}
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        atom1, atom2 = (int(value) for value in line.split()[:2])
        adjacency.setdefault(atom1, set()).add(atom2)
        adjacency.setdefault(atom2, set()).add(atom1)
    return adjacency


def _canonical_dihedral_records(text: str, adjacency: dict[int, set[int]]):
    lines = text.splitlines()
    count = int(lines[0].strip())
    records = []
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split()
        atoms = tuple(int(value) for value in fields[:4])
        tail = tuple(fields[4:])
        centers = [
            atom
            for atom in atoms
            if sum(1 for bonded in adjacency.get(atom, ()) if bonded in atoms) == 3
        ]
        if len(centers) == 1:
            center = centers[0]
            peripherals = tuple(sorted(atom for atom in atoms if atom != center))
            records.append(("periodic_improper", center, peripherals, tail))
        else:
            records.append(("exact", tuple(fields)))
    return count, Counter(records)


def _assert_spg_output_matches(actual_dir: Path, expected_dir: Path):
    actual = _snapshot_files(actual_dir)
    expected = _snapshot_files(expected_dir)

    assert actual.keys() == expected.keys()
    for name, expected_text in expected.items():
        actual_lines = actual[name].splitlines()
        expected_lines = expected_text.splitlines()
        if name == "input_dihedral.txt":
            continue
        if name == "input_angle.txt":
            assert actual_lines[0] == expected_lines[0]
            assert Counter(actual_lines[1:]) == Counter(expected_lines[1:])
            continue
        if name == "input_exclude.txt":
            assert [line.rstrip() for line in actual_lines] == [line.rstrip() for line in expected_lines]
            continue
        if name == "input.pdb":
            filtered_actual = [line for line in actual_lines if not line.startswith("CONECT")]
            filtered_expected = [line for line in expected_lines if not line.startswith("CONECT")]
            assert len(filtered_actual) == len(filtered_expected)
            for actual_line, expected_line in zip(filtered_actual, filtered_expected):
                if actual_line.startswith("CRYST1") and expected_line.startswith("CRYST1"):
                    actual_fields = actual_line.split()
                    expected_fields = expected_line.split()
                    assert actual_fields[0] == expected_fields[0] == "CRYST1"
                    for actual_field, expected_field in zip(actual_fields[1:7], expected_fields[1:7]):
                        assert abs(float(actual_field) - float(expected_field)) < 1e-6
                    assert actual_fields[7:] == expected_fields[7:]
                elif (
                    actual_line.startswith(("ATOM  ", "HETATM"))
                    and expected_line.startswith(("ATOM  ", "HETATM"))
                ):
                    assert actual_line[12:16] == expected_line[12:16]
                    assert actual_line[17:20] == expected_line[17:20]
                    assert actual_line[21] == expected_line[21]
                    for actual_field, expected_field in zip(
                        actual_line[30:54].split(),
                        expected_line[30:54].split(),
                    ):
                        assert abs(float(actual_field) - float(expected_field)) < 1e-6
                    assert actual_line[76:78].strip() == expected_line[76:78].strip()
                else:
                    assert actual_line == expected_line
            continue
        assert actual_lines == expected_lines

    adjacency = _bond_adjacency_from_text(expected["input_bond.txt"])
    actual_count, actual_records = _canonical_dihedral_records(actual["input_dihedral.txt"], adjacency)
    expected_count, expected_records = _canonical_dihedral_records(expected["input_dihedral.txt"], adjacency)
    assert actual_count == expected_count
    assert actual_records == expected_records


def _rewrite_legacy_paths(text: str, output_dir: Path) -> str:
    rewritten = text.replace(f"{LEGACY_8RYK_ROOT}/sponge", str(output_dir))
    return rewritten.replace(LEGACY_8RYK_ROOT, str(DATA_8RYK_DIR))


def _rewrite_payload_paths(value, output_dir: Path):
    if isinstance(value, str):
        return _rewrite_legacy_paths(value, output_dir)
    if isinstance(value, list):
        return [_rewrite_payload_paths(item, output_dir) for item in value]
    if isinstance(value, dict):
        return {key: _rewrite_payload_paths(item, output_dir) for key, item in value.items()}
    return value


def _run_legacy_script(script_path: Path, output_dir: Path, *, payload_override: dict | None = None,
                       include_frcmod_support: bool = False):
    env = os.environ.copy()
    if payload_override is not None:
        env["SPONGE_INTERACTIVE_FRCMOD_PAYLOAD"] = json.dumps(payload_override)
    child_code = """
import pathlib
import sys

import XpongeCPP as Xponge
import XpongeCPP.forcefield as forcefield
import XpongeCPP.forcefield.amber as amber
import XpongeCPP.forcefield.amber.ff14sb as amber_ff14sb
import XpongeCPP.forcefield.amber.gaff as amber_gaff
import XpongeCPP.forcefield.amber.tip3p as amber_tip3p

sys.modules['Xponge'] = Xponge
sys.modules['Xponge.forcefield'] = forcefield
sys.modules['Xponge.forcefield.amber'] = amber
sys.modules['Xponge.forcefield.amber.ff14sb'] = amber_ff14sb
sys.modules['Xponge.forcefield.amber.gaff'] = amber_gaff
sys.modules['Xponge.forcefield.amber.tip3p'] = amber_tip3p
if {include_frcmod_support}:
    import XpongeCPP._core as XpongeLib
    sys.modules['XpongeLib'] = XpongeLib

script_path = pathlib.Path({script_path!r})
legacy_root = {legacy_root!r}
data_root = {data_root!r}
output_dir = pathlib.Path({output_dir!r})
script_text = script_path.read_text()
script_text = script_text.replace(f"{legacy_root}/sponge", str(output_dir))
script_text = script_text.replace(legacy_root, data_root)
globals_dict = {{"__file__": str(script_path), "__name__": "__main__"}}
exec(compile(script_text, str(script_path), "exec"), globals_dict, globals_dict)
    """.format(
        include_frcmod_support="True" if include_frcmod_support else "False",
        script_path=str(script_path),
        legacy_root=LEGACY_8RYK_ROOT,
        data_root=str(DATA_8RYK_DIR),
        output_dir=str(output_dir),
    )
    result = subprocess.run(
        [sys.executable, "-c", child_code],
        check=False,
        text=True,
        capture_output=True,
        env=env,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "legacy script failed\nSTDOUT:\n"
            + result.stdout
            + "\nSTDERR:\n"
            + result.stderr
        )


def test_8ryk_spg_init_legacy_script_matches_vendored_baseline(tmp_path):
    output_dir = tmp_path / "sponge"
    output_dir.mkdir()

    _run_legacy_script(DATA_8RYK_DIR / "spg_init.txt", output_dir)

    _assert_spg_output_matches(output_dir, DATA_8RYK_DIR / "sponge")


@pytest.mark.xfail(
    reason=(
        "8RYK legacy frcmod workflow still depends on old Xponge forcefield.base "
        "and parmchk2 compatibility that XpongeCPP has not migrated yet."
    )
)
def test_8ryk_frcmod_init_legacy_script_matches_vendored_baseline(tmp_path):
    output_dir = tmp_path / "frcmod"
    output_dir.mkdir()
    payload = json.loads((DATA_8RYK_DIR / "frcmod" / "interactive.payload.json").read_text())
    payload = _rewrite_payload_paths(payload, output_dir)
    payload["boundaryMetadataPath"] = str(output_dir / "interactive.boundary.json")
    payload["rawMol2Path"] = str(output_dir / "interactive.raw.mol2")
    payload["rawOutputFrcmod"] = str(output_dir / "interactive.raw.frcmod")

    _run_legacy_script(
        DATA_8RYK_DIR / "frcmod_init.txt",
        output_dir,
        payload_override=payload,
        include_frcmod_support=True,
    )

    assert (output_dir / "interactive.raw.frcmod").read_text() == (
        DATA_8RYK_DIR / "frcmod" / "interactive.raw.frcmod"
    ).read_text()
