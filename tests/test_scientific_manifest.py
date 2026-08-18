from __future__ import annotations

from copy import deepcopy
import subprocess
import sys

from XpongeCPP.scientific_manifest import (
    _canonicalize_atom_order,
    build_scientific_manifest,
    compare_scientific_manifests,
    write_scientific_manifest,
)


def _write_case(root):
    subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import XpongeCPP as Xponge;"
                "import XpongeCPP.forcefield.amber.ff14sb;"
                "molecule=Xponge.get_peptide_from_sequence('AA');"
                "molecule.set_box_padding(4.0);"
                "Xponge.save_sponge_input("
                f"molecule,'system',{str(root)!r},format='raw')"
            ),
        ],
        check=True,
    )


def test_scientific_manifest_round_trip_and_numeric_diff(tmp_path):
    case_root = tmp_path / "case"
    _write_case(case_root)
    manifest = build_scientific_manifest(case_root, case_id="ala2")

    assert manifest["schema"] == "xponge.scientific_manifest"
    assert manifest["case_id"] == "ala2"
    assert "/atoms/mass" in manifest["datasets"]
    assert "/atoms/charge" in manifest["datasets"]
    assert "/forcefield/bond/atoms" in manifest["datasets"]
    assert "/restart/position" in manifest["datasets"]
    assert compare_scientific_manifests(manifest, manifest)["ok"]

    changed = deepcopy(manifest)
    changed["datasets"]["/atoms/charge"]["values"][0] += 1.0e-4
    diff = compare_scientific_manifests(manifest, changed)
    assert not diff["ok"]
    assert diff["issues"][0]["path"] == "/atoms/charge"

    output = write_scientific_manifest(manifest, tmp_path / "manifest.json")
    assert output.is_file()


def test_scientific_manifest_canonicalizes_improper_permutations():
    from XpongeCPP.scientific_manifest import _canonicalize_dihedrals

    def manifest_atoms(values):
        return {
            "/forcefield/bond/atoms": {
                "dtype": "int32",
                "shape": [3, 2],
                "values": [[2, 1], [2, 3], [2, 4]],
            },
            "/forcefield/dihedral/atoms": {
                "dtype": "int32",
                "shape": [1, 4],
                "values": [values],
            },
            "/forcefield/dihedral/periodicity": {
                "dtype": "int32",
                "shape": [1],
                "values": [2],
            },
            "/forcefield/dihedral/k": {
                "dtype": "float32",
                "shape": [1],
                "values": [10.5],
            },
            "/forcefield/dihedral/phi0": {
                "dtype": "float32",
                "shape": [1],
                "values": [3.141593],
            },
        }

    first = manifest_atoms([1, 3, 2, 4])
    second = manifest_atoms([4, 2, 1, 3])
    _canonicalize_dihedrals(first)
    _canonicalize_dihedrals(second)
    assert first["/forcefield/dihedral/atoms"] == second["/forcefield/dihedral/atoms"]
    assert first["/forcefield/dihedral/kind"]["values"] == ["improper"]


def test_scientific_manifest_canonicalizes_atom_insertion_order():
    first = {
        "/atoms/name": {
            "dtype": "string",
            "shape": [3],
            "values": ["C", "A", "B"],
        },
        "/atoms/type_name": {
            "dtype": "string",
            "shape": [3],
            "values": ["tc", "ta", "tb"],
        },
        "/atoms/residue_index": {
            "dtype": "int32",
            "shape": [3],
            "values": [0, 0, 0],
        },
        "/atoms/charge": {
            "dtype": "float32",
            "shape": [3],
            "values": [3.0, 1.0, 2.0],
        },
        "/restart/position": {
            "dtype": "float32",
            "shape": [3, 3],
            "values": [[3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        },
        "/forcefield/bond/atoms": {
            "dtype": "int32",
            "shape": [2, 2],
            "values": [[0, 2], [2, 1]],
        },
        "/forcefield/bond/k": {
            "dtype": "float32",
            "shape": [2],
            "values": [20.0, 10.0],
        },
        "/forcefield/bond/r0": {
            "dtype": "float32",
            "shape": [2],
            "values": [2.0, 1.0],
        },
        "/topology/exclusions/offset": {
            "dtype": "int32",
            "shape": [4],
            "values": [0, 0, 1, 1],
        },
        "/topology/exclusions/list": {
            "dtype": "int32",
            "shape": [1],
            "values": [2],
        },
    }
    second = deepcopy(first)
    second["/atoms/name"]["values"] = ["A", "B", "C"]
    second["/atoms/type_name"]["values"] = ["ta", "tb", "tc"]
    second["/atoms/charge"]["values"] = [1.0, 2.0, 3.0]
    second["/restart/position"]["values"] = [
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ]
    second["/forcefield/bond/atoms"]["values"] = [[2, 1], [1, 0]]
    second["/forcefield/bond/k"]["values"] = [20.0, 10.0]
    second["/forcefield/bond/r0"]["values"] = [2.0, 1.0]
    second["/topology/exclusions/list"]["values"] = [0]

    _canonicalize_atom_order(first)
    _canonicalize_atom_order(second)
    assert first == second


def test_scientific_manifest_canonicalizes_cmap_grid_order_and_duplicates():
    from XpongeCPP.scientific_manifest import _canonicalize_cmaps

    first = {
        "/forcefield/cmap/atoms": {
            "dtype": "int32",
            "shape": [3, 5],
            "values": [[5, 6, 7, 8, 9], [0, 1, 2, 3, 4], [10, 11, 12, 13, 14]],
        },
        "/forcefield/cmap/type": {
            "dtype": "int32",
            "shape": [3],
            "values": [0, 1, 2],
        },
        "/forcefield/cmap/resolution": {
            "dtype": "int32",
            "shape": [3],
            "values": [2, 2, 2],
        },
        "/forcefield/cmap/grid_value": {
            "dtype": "float32",
            "shape": [12],
            "values": [4.0, 3.0, 2.0, 1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0],
        },
    }
    second = {
        "/forcefield/cmap/atoms": deepcopy(first["/forcefield/cmap/atoms"]),
        "/forcefield/cmap/type": {
            "dtype": "int32",
            "shape": [3],
            "values": [1, 0, 1],
        },
        "/forcefield/cmap/resolution": {
            "dtype": "int32",
            "shape": [2],
            "values": [2, 2],
        },
        "/forcefield/cmap/grid_value": {
            "dtype": "float32",
            "shape": [8],
            "values": [1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0],
        },
    }

    _canonicalize_cmaps(first)
    _canonicalize_cmaps(second)

    assert first == second
