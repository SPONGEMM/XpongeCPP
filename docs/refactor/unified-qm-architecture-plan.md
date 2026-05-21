# Unified QM Architecture Plan

This document describes a unified quantum-mechanics (QM) architecture for
`XpongeCPP`. The goal is to turn the current RESP-specific `PySCF` / `Psi4`
+backend split into a reusable QM subsystem that can serve:

- RESP charge fitting
- future MCPB.py-like metal-center workflows
- future single-point, optimization, ESP, Hessian, and force-constant tasks

This plan is intentionally broader than the existing RESP migration plan in
[resp-psi4-cpp-migration-plan.md](/mnt/data8t/Software/Xponge/Xponge-CPP/docs/migration/resp-psi4-cpp-migration-plan.md).
That document focuses on delivering multi-backend RESP. This document focuses on
the longer-lived architecture that RESP and future QM features should share.

## Goals

- Define one reusable QM backend layer for `PySCF`, `Psi4`, and later engines.
- Make RESP a consumer of QM services rather than the owner of QM backend code.
- Create an architecture that can support MCPB.py-like workflows without
  duplicating backend integration logic.
- Keep `Assign.calculate_charge("RESP", ...)` compatible while moving backend
  logic out of `assign/`.
- Preserve the current default behavior:
  - `PySCF` remains the default RESP backend where available
  - `Psi4` remains an explicit alternative backend
- Keep numerical RESP logic independent from QM backend internals.

## Non-goals

- Reimplement SCF, gradients, or Hessians in C++.
- Replace all Python orchestration with C++.
- Commit to one permanent QM backend.
- Implement MCPB itself in this phase.

## Why The Current RESP-local Backend Layout Is Not Enough

Current RESP-specific backend modules:

- [resp.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp.py)
- [pyscf_backend.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/pyscf_backend.py)
- [psi4_backend.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/psi4_backend.py)
- [resp_core.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp_core.py)

This layout works for RESP, but it does not scale well once `XpongeCPP` needs
additional QM-driven workflows.

If future MCPB-like work is added directly on top of the current layout, the
project will likely repeat the same backend logic in multiple places:

- molecule construction
- charge and spin normalization
- SCF execution
- optional geometry optimization
- ESP computation
- Hessian and force-constant extraction
- backend-specific error handling

That would create feature-local backend layers instead of one shared QM layer.

## Target Architecture

The target architecture has four layers.

### 1. Compatibility layer

Responsibility:

- Preserve existing public APIs such as
  `Assign.calculate_charge("RESP", ...)`
- Translate legacy user parameters into the new QM API
- Keep user-facing error messages and Windows guidance stable

Examples:

- `src/XpongeCPP/_compat/assign.py`
- `src/XpongeCPP/assign/resp.py`

This layer should not directly know how `PySCF` or `Psi4` build molecules.

### 2. Feature layer

Responsibility:

- Implement feature-specific workflows such as RESP or future MCPB
- Ask the QM API for needed results
- Own feature-specific algorithms and constraints

Examples:

- `RESP` requests:
  - SCF wavefunction-like result
  - ESP on grid
- future `MCPB` requests:
  - geometry optimization
  - Hessian
  - atom/fragment selection helpers
  - force-constant extraction inputs

This layer should depend on backend-neutral QM requests and results.

### 3. QM API and scheduler layer

Responsibility:

- Offer one stable API for all QM-enabled features
- Select backend
- Validate requested capability
- Normalize inputs
- Normalize outputs
- Hide backend-specific imports and object models
- Surface consistent error handling

This is the key new layer missing today.

### 4. QM backend implementations

Responsibility:

- Convert backend-neutral requests into concrete `PySCF` / `Psi4` calls
- Build molecules
- Run SCF
- Compute ESP
- Optimize geometry
- Compute Hessians when supported

This layer contains the actual backend adapters.

## Proposed Python Module Layout

Suggested new package:

- `src/XpongeCPP/qm/__init__.py`
- `src/XpongeCPP/qm/api.py`
- `src/XpongeCPP/qm/models.py`
- `src/XpongeCPP/qm/capabilities.py`
- `src/XpongeCPP/qm/errors.py`
- `src/XpongeCPP/qm/scheduler.py`
- `src/XpongeCPP/qm/backends/base.py`
- `src/XpongeCPP/qm/backends/pyscf_backend.py`
- `src/XpongeCPP/qm/backends/psi4_backend.py`

Suggested feature consumers:

- `src/XpongeCPP/assign/resp.py`
- future `src/XpongeCPP/qm_workflows/mcpb.py`

Suggested C++ consumers:

- existing RESP numerical core in `cpp/assign/resp.cpp`
- possible future force-constant post-processing in `cpp/`

## Backend-neutral Data Models

Define backend-neutral models in `qm/models.py`.

### `QMMolecule`

Purpose:

- represent the QM system independent of backend

Suggested fields:

- `atom_symbols`
- `coordinates_angstrom`
- `total_charge`
- `spin`
- optional `atom_names`
- optional `formal_charges`
- optional `bonds`
- optional `metadata`

### `QMRunOptions`

Purpose:

- hold options common to SCF-like calculations

Suggested fields:

- `backend`
- `basis`
- `method`
- `reference`
- `optimize_geometry`
- `threads`
- `memory`
- `properties`

### `ESPGridRequest`

Purpose:

- request ESP values on explicit grid points

Suggested fields:

- `grid_points_bohr`
- optional `include_nuclear_term`
- optional `include_electronic_term`

### `SCFResult`

Purpose:

- hold common SCF outputs

Suggested fields:

- `backend_name`
- `total_energy`
- `converged`
- `coordinates_bohr`
- `nuclear_charges`
- optional backend-private handle
- optional timing metadata

### `ESPResult`

Purpose:

- return backend-neutral ESP outputs

Suggested fields:

- `grid_points_bohr`
- `electronic_esp_au`
- optional `total_esp_au`
- optional `nuclear_esp_au`
- optional timing metadata

### `OptimizationResult`

Purpose:

- support geometry optimization workflows

Suggested fields:

- `optimized_coordinates_angstrom`
- `converged`
- `iterations`
- optional `final_energy`

### `HessianResult`

Purpose:

- support future MCPB-like force-constant workflows

Suggested fields:

- `cartesian_hessian_au`
- `coordinates_angstrom`
- `atom_symbols`
- optional timing metadata

## Backend Interface Contract

Define a backend protocol in `qm/backends/base.py`.

Suggested capabilities:

- `build_molecule(qm_molecule)`
- `run_scf(qm_molecule, run_options) -> SCFResult`
- `compute_esp(scf_result, esp_request) -> ESPResult`
- `optimize_geometry(qm_molecule, run_options) -> OptimizationResult`
- `compute_hessian(qm_molecule, run_options) -> HessianResult`
- `capabilities() -> QMCapabilitySet`

Backends should not expose raw `PySCF` or `Psi4` objects outside the backend
package except as opaque private handles inside result objects.

## Capability Model

Define explicit capability flags in `qm/capabilities.py`.

Suggested flags:

- `supports_scf`
- `supports_esp`
- `supports_geometry_optimization`
- `supports_hessian`
- `supports_open_shell`
- `supports_point_charges`
- `supports_constraints`

Why this matters:

- RESP needs only `SCF + ESP`
- MCPB-like workflows will likely need
  `SCF + optimization + Hessian`
- not every backend needs to support every capability on day one

## QM Scheduler Responsibilities

The scheduler in `qm/scheduler.py` should:

- resolve backend name
- choose default backend when omitted
- check capability support before running
- lazily import backend packages
- apply platform-sensitive messages
- centralize optional dependency errors
- offer one entrypoint for features

Suggested entrypoints:

- `run_scf(...)`
- `compute_esp(...)`
- `optimize_geometry(...)`
- `compute_hessian(...)`
- `get_backend(name=None)`

## Default Backend Policy

The scheduler should own default selection logic.

### Short-term policy

- RESP defaults to `PySCF` where installed
- `Psi4` is an explicit opt-in backend
- Windows guidance remains:
  use `Psi4` explicitly when `PySCF` is unavailable

### Long-term policy

Do not hardcode all future QM features to one backend policy.

Examples:

- RESP may default to `PySCF`
- MCPB may later prefer a backend with stronger Hessian support

Therefore, default selection should allow:

- feature-specific defaults
- environment-aware defaults
- user override

## How RESP Should Use The Unified QM Layer

RESP should move from:

- importing `assign.pyscf_backend`
- importing `assign.psi4_backend`

to:

- building a backend-neutral `QMMolecule`
- requesting `SCFResult`
- requesting `ESPResult` on explicit grids
- passing backend-neutral arrays into RESP core

Long-term RESP flow:

1. adapt `Assign` -> `QMMolecule`
2. select backend through scheduler
3. generate MK grid in RESP core
4. request ESP on that grid from QM API
5. run RESP numerical fit in C++
6. write charges back to `Assign`

In this shape, RESP becomes a QM client rather than a QM backend host.

## How Future MCPB-like Work Should Use The Unified QM Layer

Future MCPB-like workflows will likely need more than RESP.

Potential MCPB-style needs:

- fragment or model-system construction
- geometry optimization
- Hessian computation
- bond / angle / force-constant extraction
- metal-center atom selection and mapping
- charge and topology integration

The important architectural point is:

MCPB should reuse the same QM API as RESP, not create a second feature-local
backend layer.

Suggested MCPB flow:

1. build selected model system from `Molecule` / `Assign`
2. adapt it to `QMMolecule`
3. request optimization and Hessian through scheduler
4. run force-constant fitting and parameter derivation in feature-specific code
5. write results into force-field artifacts

## Proposed Refactor Path

This refactor should be incremental and low-risk.

### Phase 0. Freeze current RESP state

Already available:

- multi-backend RESP routing
- C++ RESP core
- benchmark and parity tests

Use current RESP behavior as the baseline.

### Phase 1. Create `qm/` package without changing behavior

Tasks:

- add `qm/models.py`
- add `qm/errors.py`
- add `qm/capabilities.py`
- add `qm/backends/base.py`
- move backend-specific helper code conceptually into `qm/backends/`

Acceptance:

- no public behavior change
- RESP still runs through existing entrypoint

### Phase 2. Move backend code out of `assign/`

Tasks:

- create `qm/backends/pyscf_backend.py`
- create `qm/backends/psi4_backend.py`
- keep thin compatibility shims in `assign/` during transition if needed

Acceptance:

- RESP behavior remains unchanged
- backend code no longer lives primarily under `assign/`

### Phase 3. Introduce scheduler API

Tasks:

- add `qm/scheduler.py`
- add feature-neutral backend selection
- centralize import errors and platform hints

Acceptance:

- RESP obtains backend through scheduler rather than directly importing backend modules

### Phase 4. Adapt RESP to QM API

Tasks:

- `resp.py` builds `QMMolecule`
- `resp.py` requests SCF + ESP through scheduler
- `resp.py` passes neutral arrays into RESP core

Acceptance:

- no loss of current RESP features
- current RESP tests keep passing

### Phase 5. Add non-RESP QM smoke workflows

Tasks:

- create minimal API coverage for:
  - SCF only
  - optimization
  - Hessian capability detection

Acceptance:

- architecture is proven to support more than RESP

### Phase 6. Start MCPB-oriented prototype

Tasks:

- define model-system builder
- define Hessian consumer API
- prototype one small bonded-metal derivation path

Acceptance:

- first MCPB-style experiment uses the shared QM layer rather than ad hoc backend code

## Proposed Initial File Moves

Recommended first concrete move:

- move
  - `src/XpongeCPP/assign/pyscf_backend.py`
  - `src/XpongeCPP/assign/psi4_backend.py`
- to
  - `src/XpongeCPP/qm/backends/pyscf_backend.py`
  - `src/XpongeCPP/qm/backends/psi4_backend.py`

Then update:

- `src/XpongeCPP/assign/resp.py`

to import through the new `qm` layer rather than owning backend selection.

## Testing Strategy

Testing should be split by layer.

### Compatibility tests

Examples:

- `Assign.calculate_charge("RESP", ...)` still works
- default backend behavior remains stable
- Windows-facing `Psi4` guidance remains correct

### QM backend tests

Examples:

- `PySCF` can run SCF and ESP on fixture systems
- `Psi4` can run SCF and ESP on fixture systems
- geometry optimization smoke tests where available
- Hessian capability detection tests

### Feature tests

Examples:

- RESP parity
- RESP backend comparison
- RESP Python/C++ core comparison

### Architecture tests

Examples:

- scheduler rejects unsupported capability requests cleanly
- feature code can switch backend without backend-specific imports

## Risks

### Risk 1. Over-generalizing too early

If the QM API is made too abstract before real consumers exist, it may become
hard to use.

Mitigation:

- design around RESP now
- reserve space for MCPB requirements
- add only the next real capabilities needed

### Risk 2. Feature code leaking backend specifics back upward

If RESP or MCPB keeps reaching into `psi4` or `pyscf` private objects, the new
layer loses value.

Mitigation:

- keep backend-private objects opaque
- exchange plain arrays and simple result objects at boundaries

### Risk 3. Different feature defaults fighting each other

RESP and MCPB may eventually want different preferred backends.

Mitigation:

- keep defaults at the scheduler + feature-policy level
- avoid one hardcoded global backend rule

### Risk 4. Test coverage staying RESP-only

If all tests remain RESP-based, the new QM layer may still be implicitly
RESP-shaped.

Mitigation:

- add non-RESP API smoke tests early

## Immediate Recommendation

The next refactor after the current RESP migration should be:

1. create `src/XpongeCPP/qm/`
2. move `PySCF` and `Psi4` backend code under `qm/backends/`
3. add backend-neutral data models and scheduler
4. rewire `RESP` to consume the QM scheduler

That gives `XpongeCPP` a real QM subsystem while keeping RESP stable and making
future MCPB-like work much cleaner.
