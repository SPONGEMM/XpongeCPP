#!/usr/bin/env python3
"""Generate the original-Xponge 1KV2+B96 Amber reference outputs.

This script is intentionally kept outside pytest because it depends on the
local original Xponge checkout and can take seconds to minutes on large boxes.
It records the exact Python calls used as the parity reference for XpongeCPP.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


SPONGE_KEYS = [
    "residue",
    "resname",
    "atom_name",
    "atom_type_name",
    "coordinate",
    "mass",
    "charge",
    "LJ",
    "bond",
    "angle",
    "dihedral",
    "exclude",
    "nb14",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xponge-repo",
        type=Path,
        default=Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge"),
        help="Path to the original Xponge checkout.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data"),
        help="Directory containing 1KV2_H.pdb, B96.mol2, B96_H.mol2, and B96.frcmod.",
    )
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory for generated reference files.")
    parser.add_argument("--prefix", default="spg", help="SPONGE input file prefix.")
    parser.add_argument("--padding", type=float, default=10.0, help="Add_Solvent_Box distance in Angstrom.")
    parser.add_argument("--tolerance", type=float, default=2.5, help="Add_Solvent_Box solvent tolerance.")
    parser.add_argument("--seed", type=int, default=20260509, help="NumPy seed used before water and ion replacement.")
    parser.add_argument("--na", type=int, default=64, help="Number of NA residues to replace from WAT.")
    parser.add_argument("--cl", type=int, default=52, help="Number of CL residues to replace from WAT.")
    parser.add_argument(
        "--assign-b96",
        action="store_true",
        help="Use original Xponge GAFF assignment on B96_H.mol2 and save typed_b96.mol2 before assembly.",
    )
    parser.add_argument(
        "--no-solvent",
        action="store_true",
        help="Only assemble protein+B96 and export; useful for fast topology debugging.",
    )
    return parser.parse_args()


def import_original_xponge(xponge_repo: Path):
    sys.path.insert(0, str(xponge_repo))
    import Xponge  # pylint: disable=import-error,import-outside-toplevel
    import Xponge.forcefield.amber.ff14sb  # noqa: F401 pylint: disable=import-error,import-outside-toplevel
    import Xponge.forcefield.amber.gaff  # noqa: F401 pylint: disable=import-error,import-outside-toplevel
    import Xponge.forcefield.amber.tip3p  # noqa: F401 pylint: disable=import-error,import-outside-toplevel
    from Xponge import ResidueType  # pylint: disable=import-error,import-outside-toplevel
    from Xponge.forcefield import amber  # pylint: disable=import-error,import-outside-toplevel
    from Xponge.process import Add_Solvent_Box, Solvent_Replace  # pylint: disable=import-error,import-outside-toplevel

    return Xponge, ResidueType, amber, Add_Solvent_Box, Solvent_Replace


def assign_b96_with_original_xponge(Xponge, data_dir: Path, out_dir: Path) -> Path:
    source = data_dir / "B96_H.mol2"
    typed = out_dir / "typed_b96.mol2"
    assignment = Xponge.get_assignment_from_mol2(str(source), total_charge="sum")
    assignment.determine_atom_type("gaff")
    restype = assignment.to_residuetype("B")
    molecule = Xponge.Molecule(restype)
    Xponge.save_mol2(molecule, str(typed))
    return typed


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    Xponge, ResidueType, amber, Add_Solvent_Box, Solvent_Replace = import_original_xponge(args.xponge_repo)
    amber.load_parameters_from_frcmod(str(args.data_dir / "B96.frcmod"), prefix=False)

    b96_mol2 = assign_b96_with_original_xponge(Xponge, args.data_dir, args.out_dir) if args.assign_b96 else (
        args.data_dir / "B96.mol2"
    )

    molecule = Xponge.load_pdb(str(args.data_dir / "1KV2_H.pdb"))
    ligand = Xponge.load_mol2(str(b96_mol2))
    for residue in ligand.residues:
        molecule.Add_Residue(residue)

    if not args.no_solvent:
        np.random.seed(args.seed)
        Add_Solvent_Box(
            molecule,
            ResidueType.get_type("WAT"),
            args.padding,
            tolerance=args.tolerance,
        )
        np.random.seed(args.seed)
        Solvent_Replace(
            molecule,
            ResidueType.get_type("WAT"),
            {
                ResidueType.get_type("NA"): args.na,
                ResidueType.get_type("CL"): args.cl,
            },
        )

    Xponge.save_pdb(molecule, str(args.out_dir / f"{args.prefix}.pdb"))
    Xponge.save_sponge_input(molecule, prefix=args.prefix, dirname=str(args.out_dir))

    print(f"residues {len(molecule.residues)}")
    print(f"atoms {sum(len(residue.atoms) for residue in molecule.residues)}")
    for key in SPONGE_KEYS:
        path = args.out_dir / f"{args.prefix}_{key}.txt"
        print(f"{key} {path.read_text().splitlines()[0]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
