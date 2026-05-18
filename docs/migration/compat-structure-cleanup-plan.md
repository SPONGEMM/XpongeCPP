# XpongeCPP Compatibility Structure Cleanup Plan

## Goal

Reduce the structural complexity introduced by the first-wave legacy
compatibility migration, while preserving the already-green compatibility test
surface.

This cleanup focuses on architecture and maintenance cost, not on expanding
legacy feature coverage.

Legacy-script preservation rule:

- do not trade away existing old-`Xponge` script compatibility merely to make
  the package layout look cleaner
- `import Xponge` and old `Xponge.*` import paths remain first-class
  compatibility targets
- structural cleanup should prefer internal consolidation over user-visible
  behavior changes

## Current Progress

Phase 1, the first half of Phase 3, Phase 5, and the first split step of Phase
6 are now in progress or complete:

- a dedicated bootstrap module now exists at
  `src/XpongeCPP/_compat/bootstrap.py`
- automatic install side effects were removed from `_compat/api.py`
- `src/XpongeCPP/__init__.py` now calls the centralized bootstrap entrypoint
  instead of directly orchestrating several patch installers itself
- `src/XpongeCPP/compat.py` now treats `_compat.api` as the canonical public
  source of compatibility helpers
- `src/XpongeCPP/_compat/__init__.py` has been reduced to an internal module
  package index instead of re-exporting the same public helper surface again
- the highest-intrusion logic in `src/Xponge/__init__.py` has started moving
  into `_compat.imports`:
  - public attribute forwarding now uses a shared helper instead of a local
    hand-written loop
  - `__main__` namespace mirroring now uses a shared helper instead of local
    inline mutation code
- a first audited `__main__` injection policy now exists:
  - default policy remains strict legacy compatibility
  - `XPONGECPP_LEGACY_MAIN_NAMESPACE=0` disables bare-name mirroring for
    explicit opt-out/debug cases
  - existing old-script behavior remains the default and is covered by tests
- `src/XpongeCPP/mdrun.py` now writes `BIN_PATH.dat` to runtime state under the
  user cache directory instead of polluting `src/`
- `.gitignore` now ignores the old `src/XpongeCPP/BIN_PATH.dat` path and the
  stale in-tree runtime file has been removed locally
- repository-local ignore rules now also cover `src/**/__pycache__/` and
  `tests/**/__pycache__/`, and the existing in-tree bytecode caches were
  cleaned out during the current pass so the source tree is free of known
  runtime-generated artifacts
- the first three compatibility test splits have landed:
  - `tests/test_compat_surface.py` now carries the import/package-surface and
    first-wave legacy workflow skeleton checks that previously lived at the top
    of `tests/test_compat_layer.py`
  - `tests/test_compat_assign.py` now carries the old-`Assign` object-shape,
    AtomType registry, and first-wave assignment DSL compatibility checks that
    previously lived in `tests/test_compat_layer.py`
  - `tests/test_compat_process.py` now carries the executable process-compat
    checks for `AddSolventBox`, `Add_Molecule`, `AddIons`,
    `Impose_*`, `H_Mass_Repartition`, `Main_Axis_Rotate`,
    `Get_Peptide_From_Sequence`, `SetBoxPadding`, and `Solvent_Replace`
- the primary evidence command has been updated and remains green after the test
  split:
  `pixi run pytest -q tests/test_compat_surface.py tests/test_compat_assign.py tests/test_compat_process.py tests/test_compat_layer.py tests/test_legacy_import_matrix.py tests/test_8ryk_manual_frcmod.py tests/test_poly_non_amber_workflow.py tests/test_b96_mol2_gaff.py`

This is intentionally only a first cleanup step. Legacy script behavior remains
preserved; the public `Xponge` import path and existing compatibility tests
must stay green while later phases continue.

## Motivation

The current compatibility layer is functionally much stronger than before, but
it still carries several structural risks:

1. Import-time side effects are too strong.
2. Public compatibility entrypoints are partially duplicated.
3. The `src/Xponge/` shim tree is verbose and repetitive.
4. Some runtime-generated state still lands inside `src/`.
5. Compatibility tests are concentrated in a few oversized files.

These issues make the codebase harder to reason about, harder to isolate in
tests, and more expensive to evolve.

## Scope

In scope:

- `src/XpongeCPP/__init__.py`
- `src/Xponge/__init__.py`
- `src/XpongeCPP/compat.py`
- `src/XpongeCPP/_compat/*`
- `src/Xponge/*` shim package layout
- runtime output placement such as `BIN_PATH.dat`
- compatibility test organization

Out of scope:

- new legacy feature coverage
- C++ core refactors unrelated to compatibility bootstrap
- `analysis` feature expansion
- `MindSponge` feature expansion

## Problems To Fix

### 1. Import-time global side effects

Current state:

- `src/XpongeCPP/__init__.py` installs runtime patches during import
- `src/Xponge/__init__.py` injects names into `__main__` during import

Risks:

- surprising global mutation
- test order sensitivity
- harder debugging in REPL / notebooks / scripts
- less explicit activation of compatibility behavior

Desired state:

- importing `XpongeCPP` should expose APIs, not silently mutate global runtime
- importing `Xponge` should stay as lightweight as practical **without**
  breaking existing old-script expectations unless an audit proves that the
  behavior is unused
- high-intrusion namespace injection must be audited before being weakened or
  removed

### 2. Duplicated compatibility public surface

Current state:

- `src/XpongeCPP/_compat/__init__.py`
- `src/XpongeCPP/compat.py`

Both expose nearly the same public compatibility helpers.

Risks:

- drift between two public entrypoints
- unclear ownership
- unnecessary maintenance overhead

Desired state:

- one canonical public compatibility entrypoint
- `_compat/` remains clearly internal

### 3. Large shim tree under `src/Xponge/`

Current state:

- many very small files only re-export from `XpongeCPP`

Risks:

- repetitive boilerplate
- wide file surface for trivial changes
- harder review and navigation

Desired state:

- reduce boilerplate where possible
- keep package-name compatibility as a hard requirement
- centralize forwarding logic without requiring old scripts to switch from
  `import Xponge` to `import XpongeCPP as Xponge`

### 4. Runtime artifacts inside `src/`

Current state:

- `src/XpongeCPP/BIN_PATH.dat` can be generated locally

Risks:

- source tree pollution
- accidental packaging concerns
- confusing repository state

Desired state:

- runtime state moved outside `src/`
- clear local cache / runtime directory policy

### 5. Oversized compatibility tests

Current state:

- `tests/test_compat_layer.py` covers many unrelated areas

Risks:

- difficult failures
- poor locality when editing behavior
- higher maintenance cost

Desired state:

- split compatibility tests by topic
- keep evidence clear and discoverable

## Proposed Refactor

### Phase 1. Separate bootstrap from API export

Create a dedicated bootstrap module, for example:

- `src/XpongeCPP/_compat/bootstrap.py`

Move import-time installation logic out of `src/XpongeCPP/__init__.py`:

- runtime patch installation
- assign patch installation
- template/global sync bootstrap
- MindSponge todo bootstrap
- optional legacy namespace wiring

Keep `XpongeCPP.__init__` focused on:

- stable symbol exports
- thin wrappers
- explicit bootstrap call only if truly required

Acceptance:

- `XpongeCPP.__init__` becomes substantially smaller
- patch-install responsibilities are concentrated in one bootstrap module

### Phase 2. Make `Xponge` shim import less intrusive

Refactor `src/Xponge/__init__.py` so that:

- first perform an audit of real old-script dependence on default `__main__`
  injection
- only after that audit, choose one of:
  - keep default injection for strict legacy mode
  - gate it behind an explicit compatibility mode
  - or narrow it without breaking known script patterns

Compatibility note:

- because the explicit user goal is “change old Xponge scripts as little as
  possible”, this phase must optimize for script preservation first and cleanup
  second
- no cleanup step in this phase is allowed to break:
  - `import Xponge`
  - old `Xponge.*` package imports
  - bare-name script patterns already covered by compatibility tests

Acceptance:

- an audited policy exists for `__main__` injection
- existing package-name import compatibility remains intact
- no regression in known bare-name old-script patterns

### Phase 3. Canonicalize public compatibility entrypoint

Pick one official public compatibility module:

- recommended: `src/XpongeCPP/compat.py`

Then:

- make `src/XpongeCPP/_compat/__init__.py` internal-only
- avoid exporting overlapping public API from both places
- optionally keep `_compat.__init__` minimal for internal imports only

Acceptance:

- one documented public compatibility entrypoint
- no duplicated public compatibility export surface

### Phase 4. Reduce shim boilerplate in `src/Xponge/`

Audit the shim files and group them into categories:

1. simple re-export only
2. package-path bridge only
3. special-case shim with custom behavior

For category 1 and 2, prefer more centralized forwarding patterns where safe.

Possible directions:

- shared forwarding helpers
- thinner per-package `__init__` files
- generated shim tables if maintenance cost remains high

Acceptance:

- lower repetitive shim boilerplate
- no regression in `test_legacy_import_matrix.py`
- no requirement for users to rewrite old imports as
  `import XpongeCPP as Xponge`

### Phase 5. Move runtime state out of `src/`

Relocate files like `BIN_PATH.dat` to one of:

- repository runtime temp dir
- user cache dir
- configurable app-state path

Add explicit ignore policy if needed.

Acceptance:

- no runtime-generated artifacts inside `src/`
- local reruns do not dirty the source tree with runtime state

### Phase 6. Split compatibility tests by concern

Refactor the current large compatibility test into smaller files such as:

- `tests/test_compat_imports.py`
- `tests/test_compat_runtime.py`
- `tests/test_compat_assign.py`
- `tests/test_compat_process.py`
- `tests/test_compat_workflows.py`
- `tests/test_compat_forcefield_base.py`

Keep the current end-to-end evidence, but improve locality.

Acceptance:

- failures clearly identify the compatibility subsystem involved
- no loss of evidence coverage

## Execution Order

1. Phase 1: bootstrap separation
2. Phase 2: audit and possibly narrow `__main__` injection without breaking
   known old-script patterns
3. Phase 3: canonicalize compatibility public API
4. Phase 4: shrink `src/Xponge/` shim boilerplate
5. Phase 5: remove runtime artifacts from `src/`
6. Phase 6: split compatibility tests

## Validation Plan

At minimum after each phase:

```bash
pixi run pytest -q tests/test_compat_surface.py tests/test_compat_assign.py tests/test_compat_process.py tests/test_compat_layer.py tests/test_legacy_import_matrix.py
pixi run pytest -q tests/test_8ryk_manual_frcmod.py tests/test_poly_non_amber_workflow.py tests/test_b96_mol2_gaff.py
```

After test splitting, update the command set to the new test file names while
preserving equivalent coverage.

## Success Criteria

This cleanup is successful when:

- compatibility bootstrap is explicit and centralized
- `XpongeCPP.__init__` becomes thinner and easier to audit
- `Xponge` shim behavior is either preserved or narrowed only after script
  dependency audit
- only one public compatibility entrypoint remains
- runtime artifacts no longer appear in `src/`
- compatibility tests are split into maintainable topical files
- no regression in the current green legacy compatibility evidence
- old scripts still prefer and retain `import Xponge` / `Xponge.*` paths
