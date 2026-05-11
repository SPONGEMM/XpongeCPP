---
name: xpongecpp-gaff-assign-migration
description: Use when migrating, extending, or testing XpongeCPP Assign and Amber GAFF atom typing against original Xponge, especially for mol2 inputs, AssignRule parity, ring markers, bond markers, and 100-molecule validation sets.
---

# XpongeCPP GAFF Assign Migration

## Core Rule

Port original Xponge Assign/GAFF behavior. Do not add molecule-name, atom-name, or one-off chemistry shortcuts to pass a single molecule.

Reference sources:

- `/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge/Xponge/assign/__init__.py`
- `/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge/Xponge/assign/__init__.py` `_RING` and `determine_ring_and_bond_type`
- `/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge/Xponge/assign/bond_order.py`
- `/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge/Xponge/forcefield/amber/gaff.py`

## Architecture

Assign must live outside core molecule storage:

- `cpp/assign/assign.cpp`: `Assign` state and public API.
- `cpp/assign/mol2_assignment.cpp`: `get_assignment_from_mol2_text`.
- `cpp/assign/ring.cpp`: ring detection, ring size markers, aromatic markers.
- `cpp/assign/gaff.cpp`: GAFF rule registry and rule functions.
- `cpp/assign/bond_order.cpp`: Xponge bond-order assignment, added incrementally.

`cpp/core/` should keep only the flat molecular data model. Python bindings may expose Assign, but hot loops and rule evaluation stay in C++.

## Migration Order

1. Preserve current B96 parity tests before refactoring.
2. Move Assign code into `cpp/assign/` without behavior changes.
3. Rebuild and run B96/1KV2/full tests.
4. Port original marker generation:
   - `atom_judge`;
   - `bond_marker` and `atom_marker`;
   - `determine_ring_and_bond_type`;
   - ring size markers `RG`, `RG3`, `RG4`, `RG5`, `RG6`;
   - aromatic markers `AR1`, `AR2`, `AR3`, `AR4`, `AB`.
5. Port GAFF rules in original registration order.
6. Add a rule coverage test that fails if any original `@gaff.Add_Rule(...)` is missing in C++.
7. Add 100-molecule atom-type parity against original Xponge baseline.

## 100-Molecule Validation

Use a fixed manifest under `tests/data/gaff_assign_100/manifest.json`.

Each entry should include:

- source id from a fixed small-molecule database snapshot/query, preferably ChEMBL ID or PubChem CID;
- input mol2 path;
- original Xponge GAFF atom type sequence;
- optional original typed mol2 path;
- optional notes if original Xponge cannot type or parameterize it.

Tests compare atom type sequences exactly. SPONGE topology comparison is only required for molecules with complete parameters/frcmod.

## Required Tests

- B96 `B96_H.mol2 -> determine_atom_type("gaff")` stays identical to original Xponge.
- Typed mol2 roundtrip preserves atom name, atom type, charge, coordinate, residue, and connectivity.
- 100-molecule manifest test compares XpongeCPP atom types to stored original Xponge atom types.
- Rule coverage test compares original `gaff.py` rule names to C++ registered rule names.
- Full pytest must pass after every migration slice.

## Red Flags

- Any branch using a specific molecule name like `B96`.
- Any branch using specific atom names such as `C1`, `N44`, etc.
- Rule order differs from original Xponge without documented proof.
- Missing rule silently falls back to element names.
- Ring/aromatic behavior is guessed instead of derived from Xponge's ring logic.
