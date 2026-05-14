# Xponge-CPP Code Structure Reorganization Plan

## Objective

Reorganize the current code structure without changing user-visible behavior. The goal is to split oversized translation units and the top-level Python facade into clearer modules while preserving:

- `import XpongeCPP as Xponge`
- the current `8ryk` `spg_init` behavior
- current `1kv2` workflow behavior
- no intentional semantic changes
- no obvious benchmark regression

## Execution Rules

These rules are mandatory for every implementation round:

1. Before making any code change, reread this Markdown file.
2. Keep each phase as a pure refactor unless this file is explicitly amended first.
3. After each phase, update this file with:
   - completed items
   - remaining items
   - any boundary adjustments that were required
4. Run the fixed validation gate after each phase.
5. Commit each phase separately.
6. Do not move to the next phase until the current phase is green.

## Fixed Validation Gate

Run these commands after every phase:

```bash
rtk pixi run test -q -rs
rtk pixi run pytest -q tests/test_8ryk_regression.py -k spg_init
rtk pixi run python benchmarks/bench_1kv2.py --padding 8 20 --repeat 5
```

Then run targeted checks for the files touched in the current phase.

## Baseline Notes

- Branch at plan start: `codex/xpongecpp-v1`
- The plan file did not previously exist.
- `docs/refactor/` did not previously exist.
- Fixed validation gate baseline recorded on 2026-05-14.
- `rtk pixi run test -q -rs`: passed, with one expected skipped test:
  - `tests/test_process_migration.py:392` because local SPONGE executables are not available
- `rtk pixi run pytest -q tests/test_8ryk_regression.py -k spg_init`: passed
- `rtk pixi run python benchmarks/bench_1kv2.py --padding 8 20 --repeat 5`: passed
- `bench_1kv2.py` median wall time snapshot:
  - `padding=8.0A total=0.088330s`
  - `padding=20.0A total=0.183149s`
- Current large files targeted by this reorganization:
  - `cpp/core/molecule.cpp`
  - `cpp/forcefield/amber.cpp`
  - `cpp/forcefield/nonamber.cpp`
  - `cpp/io/pdb.cpp`
  - `cpp/python/bindings.cpp`
  - `src/XpongeCPP/__init__.py`

## Phase 0

Create this plan file and lock the execution boundary.

Deliverables:

- create `docs/refactor/code-structure-reorg-plan.md`
- record validation commands and phase boundaries
- run the fixed validation gate as a baseline
- commit only the plan file if baseline is green

Allowed edits:

- `docs/refactor/code-structure-reorg-plan.md`

Forbidden in this phase:

- source refactors
- public API changes
- moving implementation files

## Phase 1

Split `cpp/core/molecule.cpp` into focused implementation units.

Planned new files:

- `cpp/core/element_mass.cpp`
- `cpp/core/template_ops.cpp`
- `cpp/forcefield/special/gb.cpp`

Target moves:

- element and mass helpers into `element_mass.cpp`
- template rebuild and add-missing-atoms support into `template_ops.cpp`
- GB radius and related special-forcefield helpers into `gb.cpp`

Keep in `molecule.cpp`:

- core `Molecule` behavior
- residue and residue-type behavior
- atom-level structural editing directly tied to core data structures

Constraints:

- no solvation logic moves in this phase
- no topology logic moves in this phase
- no Python binding changes in this phase

## Phase 2

Split `cpp/forcefield/amber.cpp`.

Planned new files:

- `cpp/forcefield/amber_registry.cpp`
- `cpp/forcefield/amber_parser.cpp`
- `cpp/forcefield/amber_templates.cpp`
- `cpp/forcefield/amber_builtins.cpp`
- `cpp/forcefield/amber_names.cpp`

Constraints:

- preserve registration order
- preserve template precedence
- preserve current import side effects
- preserve public Python entrypoints

## Phase 3

Split `cpp/forcefield/nonamber.cpp`.

Planned new files:

- `cpp/forcefield/gromacs_parser.cpp`
- `cpp/forcefield/charmm_loader.cpp`
- `cpp/forcefield/opls_loader.cpp`
- `cpp/forcefield/special_pairwise.cpp`

Constraints:

- preserve loader names
- preserve supported file-format coverage
- preserve CHARMM, OPLS, SW, and EDIP behavior

## Phase 4

Split `cpp/io/pdb.cpp` and `cpp/python/bindings.cpp`.

Planned new files:

- `cpp/io/pdb_reader.cpp`
- `cpp/io/pdb_writer.cpp`
- `cpp/python/bindings_core.cpp`
- `cpp/python/bindings_io.cpp`
- `cpp/python/bindings_forcefield.cpp`
- `cpp/python/bindings_assign.cpp`

Constraints:

- keep `bindings.cpp` as the thin module aggregator
- preserve all current Python names
- preserve current PDB read and write behavior

## Phase 5

Split `src/XpongeCPP/__init__.py` and separate runtime code from data assets.

Planned new files:

- `src/XpongeCPP/process.py`
- `src/XpongeCPP/template_ops.py`
- `src/XpongeCPP/io_compat.py`
- `src/XpongeCPP/legacy_types.py`

Constraints:

- keep `src/XpongeCPP/__init__.py` as a facade and re-export surface
- preserve old import style
- preserve `8ryk` `spg_init`
- do not break current data-path usage

## Progress Log

- [x] Phase 0 plan file created
- [x] Phase 0 baseline validation complete
- [ ] Phase 0 committed
- [ ] Phase 1 complete
- [ ] Phase 2 complete
- [ ] Phase 3 complete
- [ ] Phase 4 complete
- [ ] Phase 5 complete
