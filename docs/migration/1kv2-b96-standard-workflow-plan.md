# 1KV2+B96 Standard Workflow Migration Plan

> For agentic workers: use `.codex/skills/xpongecpp-code-standards` and `.codex/skills/xpongecpp-migration-checklist` before implementation. Track steps with checkboxes. Topology behavior must come from original Xponge algorithms, not new chemistry heuristics.

**Goal:** Implement the standard Xponge-style 1KV2+B96 workflow in XpongeCPP: ff14SB protein, B96 GAFF assignment/typed mol2 path, B96 frcmod, TIP3P solvation, ion replacement, PDB export, and SPONGE input export with Xponge-equivalent topology.

**Architecture:** Keep C++ as the only source of truth. Python only exposes compatibility entry points. Use original Xponge only to generate baselines and identify algorithms.

**Tech Stack:** C++17, pybind11, CMake, scikit-build-core, pytest, original Xponge reference repository.

---

## Current Baseline

Already implemented:

- `load_pdb(1KV2_H.pdb)` for current 1KV2 protein.
- `import XpongeCPP.forcefield.amber.ff14sb`.
- `import XpongeCPP.forcefield.amber.tip3p`.
- `import XpongeCPP.forcefield.amber.gaff` and `gaff2`.
- `load_mol2(B96.mol2)` preserving explicit mol2 bonds.
- `load_frcmod(B96.frcmod)` registering Amber parameters.
- B96 single-molecule SPONGE parity against Xponge.
- 1KV2 protein + TIP3P + NA/CL workflow headers.

Known blockers before full 1KV2+B96:

- `Add_Ions` rebuilds `Molecule` and currently does not preserve `explicit_bonds`/`residue_links`.
- No system assembly API exists for protein + ligand merging.
- B96 assignment from untyped `B96_H.mol2` is not implemented in C++.
- Water/ion handling remains WAT/NA/CL-centric.
- Missing topology parameters can still be hidden by fallback defaults.

---

## Task 1: Lock the Xponge Reference Workflow

**Files:**

- Create: `tests/reference/test_xponge_1kv2_b96_reference.py` or a local script under `benchmarks/`.
- Modify only tests/docs first.

- [ ] Write a reference script using original Xponge:

```python
from pathlib import Path
import numpy as np
import Xponge
import Xponge.forcefield.amber.ff14sb
import Xponge.forcefield.amber.gaff
import Xponge.forcefield.amber.tip3p
from Xponge.forcefield import amber

base = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data")
out = Path("/tmp/xponge_1kv2_b96/original")
out.mkdir(parents=True, exist_ok=True)

np.random.seed(20260509)
amber.load_parameters_from_frcmod(str(base / "B96.frcmod"), prefix=False)
protein = Xponge.load_pdb(str(base / "1KV2_H.pdb"))
ligand = Xponge.load_mol2(str(base / "B96.mol2"))

# Use original Xponge's actual assembly operation here after identifying it.
# Then call Add_Solvent_Box, Add_Ions, save_sponge_input, save_pdb.
```

- [ ] Identify the exact Xponge API used to merge/add ligand into protein.
- [ ] Identify the exact Xponge workflow for B96 assignment from `B96_H.mol2`, including output typed mol2 generation.
- [ ] Save all 13 SPONGE files and assembled PDB.
- [ ] Record headers/counts in a test fixture.

Acceptance:

- Reference run produces deterministic outputs under `/tmp/xponge_1kv2_b96/original`.
- The script uses original Xponge APIs, not manual reconstruction.

---

## Task 2: Preserve Explicit Bonds Through Ion Replacement

**Files:**

- Modify: `cpp/solvation/solvation.cpp`
- Test: `tests/test_b96_mol2_gaff.py` or new `tests/test_system_assembly.py`

- [ ] Write a failing test:
  - load B96 mol2 + frcmod;
  - add several WAT residues;
  - call `Add_Ions`;
  - export B96 bond/topology;
  - assert B96 explicit mol2 bond count remains unchanged.

- [ ] Implement ID remapping during `Add_Ions` rebuild:
  - create `old_atom_id -> new_atom_id`;
  - copy all explicit bonds whose atoms still exist;
  - copy residue links whose atoms still exist;
  - skip bonds inside removed waters;
  - add no explicit bonds for one-atom ions.

- [ ] Run targeted test.
- [ ] Run full pytest.

Acceptance:

- `Molecule::validate()` passes after ion replacement.
- B96 explicit connectivity survives water-to-ion replacement.

---

## Task 3: Add a System Assembly API

**Files:**

- Modify: `cpp/core.hpp`
- Modify: `cpp/core/molecule.cpp`
- Modify: `cpp/python/bindings.cpp`
- Modify: `src/XpongeCPP/__init__.py`
- Test: `tests/test_system_assembly.py`

- [ ] Add `append_molecule(Molecule& target, const Molecule& source)`.
- [ ] Preserve atom order, residue order, coordinates, box state policy, explicit bonds, and residue links.
- [ ] Add Python aliases such as `Add_Molecule` only if they match a common Xponge usage.
- [ ] Test protein + B96 assembly:
  - protein atoms first;
  - ligand atoms after protein;
  - B96 residue retained;
  - B96 explicit bonds shifted correctly.

Acceptance:

- Assembled molecule validates.
- B96 single-molecule topology embedded in assembled molecule matches standalone B96 topology on shifted atom indices.

---

## Task 4: Strict Parameter Mode for Parity

**Files:**

- Modify: `cpp/topology/topology.cpp`
- Possibly modify: `cpp/core.hpp`
- Test: `tests/test_safety_edges.py`

- [ ] Add a strict parameter path for SPONGE export and parity tests.
- [ ] Replace fallback bond `300.0/distance` with an error in strict mode.
- [ ] Replace fallback angle `40.0/1.91` with an error in strict mode.
- [ ] Ensure current tests load required force fields and pass strict mode.

Acceptance:

- Missing bond/angle parameters fail loudly.
- No Xponge parity test can pass by silently guessing parameters.

---

## Task 5: C++ GAFF Assign for B96_H.mol2

**Files:**

- Create: `cpp/assign/gaff.cpp` if splitting assign from `cpp/core/molecule.cpp`.
- Modify: `cpp/core.hpp`
- Modify: `CMakeLists.txt`
- Modify: `cpp/python/bindings.cpp`
- Test: `tests/test_b96_assign.py`

- [ ] Use original Xponge/Amber GAFF assignment as reference.
- [ ] Read Xponge assign code completely before implementation.
- [ ] Support this first slice:
  - mol2 atom graph;
  - bond order from mol2;
  - aromatic markers;
  - atom type assignment for all B96_H atoms;
  - charge preservation or assignment behavior matching baseline.
- [ ] Generate typed mol2 output from C++ if needed.
- [ ] Compare C++ assigned B96 typed mol2 with Xponge-generated typed mol2:
  - atom count;
  - atom names;
  - atom types;
  - charges if assigned by workflow;
  - bonds.

Acceptance:

- `B96_H.mol2 -> Assign -> typed B96` matches original Xponge for B96.
- No hand-coded B96 atom index table is allowed.

---

## Task 6: Full 1KV2+B96 Workflow

**Files:**

- Test: `tests/test_1kv2_b96_workflow.py`
- Possibly modify: API and implementation files from previous tasks.

- [ ] Implement XpongeCPP workflow:

```python
import XpongeCPP as Xponge
import XpongeCPP.forcefield.amber.ff14sb
import XpongeCPP.forcefield.amber.gaff
import XpongeCPP.forcefield.amber.tip3p

base = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data")
Xponge.load_frcmod(str(base / "B96.frcmod"))
protein = Xponge.load_pdb(str(base / "1KV2_H.pdb"))
ligand = Xponge.load_mol2(str(base / "B96.mol2"))
mol = Xponge.Merge_Molecules([protein, ligand])
water = Xponge.get_template_molecule("WAT")
Xponge.Add_Solvent_Box(mol, water, 10.0, tolerance=2.5, seed=20260509)
Xponge.Add_Ions(mol, {"NA": 64, "CL": 52}, seed=20260509)
Xponge.Save_SPONGE_Input(mol, prefix="spg", dirname="/tmp/xpongecpp_1kv2_b96")
Xponge.save_pdb(mol, "/tmp/xpongecpp_1kv2_b96/spg.pdb")
```

- [ ] Compare against Xponge reference:
  - headers;
  - atom/residue counts;
  - B96 atom types and charges;
  - bond/angle/dihedral/exclude/nb14 numeric records.

Acceptance:

- All core files are numerically equivalent under documented canonicalization.
- Output PDB reloads in XpongeCPP with same atom/residue counts.

---

## Task 7: Generalize Solvent/Ion Support

**Files:**

- Modify: `src/XpongeCPP/forcefield/amber/`
- Modify: `cpp/solvation/solvation.cpp`
- Test: `tests/test_custom_solvent.py`

- [ ] Add public imports for `spce`, `tip4p`, `tip4pew`, `opc` only after their topology semantics are supported.
- [ ] Add virtual site/EP support before claiming TIP4P/OPC parity.
- [ ] Replace WAT-only ion replacement with selected solvent residue replacement.
- [ ] Support user custom solvent from mol2 + frcmod.

Acceptance:

- Custom one-residue solvent copies atom names/types/charges/coordinates/connectivity.
- TIP3P remains unchanged.
- TIP4P/OPC are not exposed as supported until virtual-site output is correct.

---

## Task 8: Performance Benchmark

**Files:**

- Create/modify: `benchmarks/bench_1kv2_b96.py`

- [ ] Benchmark original Xponge vs XpongeCPP:
  - PDB load;
  - mol2 load/assign;
  - assembly;
  - solvation;
  - ion replacement;
  - topology build/export.
- [ ] Repeat each benchmark at least 5 times.
- [ ] Report median wall time.

Acceptance:

- C++ hot paths are measured.
- Any path under 10x has a concrete profiling target.

---

## Final Verification

Run:

```bash
rtk .venv/bin/python -m pytest -q
```

Optional ASan/UBSan build:

```bash
rtk cmake -S . -B build/asan -DXPONGECPP_SANITIZE=ON
rtk cmake --build build/asan
```

Commit after each completed task with a focused message.
