# RESP Multi-backend and C++ Migration Plan

This plan keeps the current `PySCF` RESP path as the default backend, adds a
`Psi4` backend for cross-platform RESP workflows, and moves the RESP fitting
algorithm into `XpongeCPP` C++ over time, while keeping Python as a thin
compatibility surface for `Assign.calculate_charge("RESP", ...)`.

## Goals

- Keep `PySCF` available as the default RESP backend.
- Add `Psi4` as an alternative RESP backend, especially for Windows users.
- Move the RESP fitting algorithm from Python into `XpongeCPP` C++ for speed.
- Keep the public Python API stable for current `Assign.calculate_charge(...)`
  callers.
- Make backend selection explicit and predictable from Python.
- Validate migration on a real Mokda fixture:
  [CRO_1.charge-preview.capped.mol2](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/data/mokda_resp/CRO_1.charge-preview.capped.mol2)

## Non-goals

- Reimplement SCF or geometry optimization in C++.
- Change the existing `TPACM4` or `Gasteiger` charge paths.

## Current State

Current RESP logic lives in:

- [resp.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp.py)
- [_compat/assign.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/assign.py)
- [test_assign_charge_models.py](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/test_assign_charge_models.py)

The current implementation mixes three concerns in one Python module:

1. QM backend setup and SCF execution via `PySCF`
2. MK grid generation and ESP matrix construction
3. RESP stage-1 / stage-2 restrained fitting

That coupling makes it hard to:

- support multiple QM backends cleanly
- optimize performance
- test the numerical RESP core independently from the QM backend

## Target Architecture

Split RESP into three layers:

### 1. Python compatibility layer

Responsibility:

- keep `Assign.calculate_charge("RESP", ...)`
- normalize user parameters
- normalize backend selection such as `backend="pyscf"` or `backend="psi4"`
- import optional backend dependencies lazily
- pass plain arrays and metadata into the backend/core
- convert backend/core failures into user-facing exceptions
- issue a clear `Psi4` installation hint on Windows when `PySCF` is not a usable option

Suggested location:

- `src/XpongeCPP/_compat/assign.py`
- `src/XpongeCPP/assign/resp.py` as a thin orchestration wrapper

### 2. QM backend layers

Responsibility:

- build a QM molecule from `Assign`
- run `RHF` or `UHF`
- optionally optimize geometry
- compute ESP on a user-provided grid
- return backend-neutral data buffers

Suggested modules:

- `src/XpongeCPP/assign/pyscf_backend.py`
- `src/XpongeCPP/assign/psi4_backend.py`

Suggested backend-neutral return payload:

- `atom_symbols`
- `atom_coordinates_bohr`
- `nuclear_charges`
- `grid_points_bohr`
- `esp_values_au`
- `total_charge`
- `spin`
- optional optimized coordinates

### 3. C++ RESP core

Responsibility:

- MK grid generation
- matrix assembly for ESP fitting
- stage-1 RESP restrained solve
- stage-2 grouped solve
- equivalence constraint application

Suggested implementation areas:

- `cpp/assign/resp.cpp`
- `cpp/python/` binding entry for RESP core

Suggested pybind entrypoints:

- `generate_resp_mk_grid(...)`
- `fit_resp_from_esp(...)`
- or a single high-level:
  `run_resp_core(...)`

## Functional Migration Boundary

The clean ownership boundary should be:

- the selected QM backend computes molecular ESP over a supplied grid
- `XpongeCPP` computes grid generation and RESP fitting

This means the long-term Python flow should be:

1. build `Assign`
2. generate MK grid in C++
3. ask the selected backend for ESP values on that grid
4. run RESP fitting in C++
5. write charges back onto `Assign`

This boundary keeps backend-specific logic external to the `XpongeCPP`
numerical core.

## Function-by-function Migration Map

Current Python helpers in [resp.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp.py):

### Move to backend modules

- `_get_pyscf_mol`
- `fun.make_rdm1()` / `df.incore.aux_e2(...)`-style ESP evaluation logic

Replacement shape:

- `_build_pyscf_molecule(...)`
- `_run_pyscf_scf(...)`
- `_compute_pyscf_esp_on_grid(...)`
- `_build_psi4_molecule(...)`
- `_run_psi4_scf(...)`
- `_compute_psi4_esp_on_grid(...)`

### Move to C++

- `_fibonacci_grid`
- `_get_mk_grid`
- `_force_equivalence_q`
- `_resp_scf_kernel`
- `_find_tofit_second`
- `_atom_judge`
- `_correct_extra_equivalence`
- `_get_a20_and_b20`
- stage-1 matrix setup and solve in `resp_fit`
- stage-2 grouped solve in `resp_fit`

### Keep in Python

- lazy optional dependency import and error text
- parameter defaults and validation
- `Assign` object adaptation
- writing optimized coordinates back to `Assign` if geometry optimization is enabled

## Proposed File Layout

### Python

- `src/XpongeCPP/assign/resp.py`
  - orchestration only
- `src/XpongeCPP/assign/pyscf_backend.py`
  - default `PySCF` QM layer
- `src/XpongeCPP/assign/psi4_backend.py`
  - `Psi4`-specific QM layer

### C++

- `cpp/assign/resp.cpp`
  - RESP numerical core
- `cpp/assign/resp.hpp`
  - structs and declarations if needed
- `cpp/python/...`
  - pybind exposure

### Tests and fixtures

- [tests/test_assign_charge_models.py](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/test_assign_charge_models.py)
- [tests/data/mokda_resp/CRO_1.charge-preview.capped.mol2](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/data/mokda_resp/CRO_1.charge-preview.capped.mol2)
- [tests/data/mokda_resp/formamide_resp.mol2](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/data/mokda_resp/formamide_resp.mol2)

## Execution Phases

### Phase 0. Freeze baseline

Deliverables:

- document this plan
- ensure the Mokda preview fixture is checked in
- keep the current `PySCF` implementation intact while building comparison tests
- define the intended default backend behavior

Acceptance:

- plan file exists
- baseline RESP tests still exist for current Python implementation

### Phase 1. Add backend-independent RESP tests

Goal:

Create tests that describe RESP behavior without baking in `PySCF` internals.

Add tests for:

- total charge conservation
- result vector length equals atom count
- expected sign pattern on known atoms
- stable results for `formamide_resp.mol2`
- stable results for `CRO_1.charge-preview.capped.mol2`

Add explicit fixture-driven comparison helpers for:

- `PySCF` vs `Psi4`
- Python RESP core vs C++ RESP core
- backend selection and default behavior
- Windows-facing backend error messaging

Acceptance:

- tests can compare two RESP implementations using the same input `Assign`
- no production code has switched backend yet

### Phase 2. Introduce Psi4 backend in Python

Goal:

Add `Psi4` as an alternative QM backend while keeping `PySCF` as the default
QM backend and keeping the RESP fitting code in Python for initial parity.

Tasks:

- add `psi4_backend.py`
- implement molecule build, SCF, optional optimization, ESP-on-grid
- route RESP by explicit backend selection
- keep `PySCF` as the default backend
- on Windows, provide a clear message that `Psi4` is the supported fallback when
  `PySCF` is unavailable or unsupported
- preserve current public Python signature

Acceptance:

- `PySCF` remains the default backend
- `Psi4` path reproduces `PySCF` RESP charges within agreed tolerance
- Mokda preview fixture passes under `Psi4`

### Phase 3. Move RESP numerical core to C++

Goal:

Keep both QM backends available while migrating RESP numerical work to C++.

Tasks:

- implement MK grid generation in C++
- implement stage-1 solve in C++
- implement stage-2 grouping and solve in C++
- expose pybind API
- make Python orchestration call C++ core for both backend payloads

Acceptance:

- C++ RESP core reproduces Python RESP core within agreed tolerance
- Mokda preview fixture passes with C++ core enabled

### Phase 4. Finalize multi-backend RESP behavior

Goal:

Finalize the long-term user-facing RESP backend contract.

Tasks:

- document `PySCF` as the default backend
- document `Psi4` as the alternative backend
- document Windows-specific guidance for RESP
- ensure benchmark and comparison coverage remains in place

Acceptance:

- RESP defaults to `PySCF`
- RESP can be switched to `Psi4`
- Windows guidance is documented and tested at the error-message level

## Test and Benchmark Strategy

Primary real-world regression fixture:

- [CRO_1.charge-preview.capped.mol2](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/data/mokda_resp/CRO_1.charge-preview.capped.mol2)

Secondary small-molecule fixture:

- [formamide_resp.mol2](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/data/mokda_resp/formamide_resp.mol2)

### A. Result comparison: Psi4 vs PySCF

Purpose:

- validate backend parity before relying on `Psi4` as the Windows-compatible path

Test setup:

- same input `Assign`
- same basis, charge, spin
- same grid density and grid layer
- same `only_esp` / `two_stage` switches

Record:

- full charge vector
- total charge
- max absolute per-atom delta
- RMS delta over all atoms

Recommended tolerances for first cut:

- total charge difference <= `1e-5`
- max absolute per-atom delta <= `5e-4`
- RMS delta <= `1e-4`

If these are too strict in practice, loosen only after capturing evidence on the
Mokda and formamide fixtures.

### B. Speed comparison: Psi4 vs PySCF

Purpose:

- measure QM backend runtime impact

Benchmark cases:

- `formamide_resp.mol2`
- `CRO_1.charge-preview.capped.mol2`

Measure separately:

- molecule build time
- SCF time
- ESP-on-grid time
- total wall time

Benchmark rules:

- fixed basis, charge, spin
- fixed grid parameters
- same machine, same environment
- at least 3 repetitions after one warm-up run

Primary reporting metrics:

- median wall time
- min/max wall time
- speed ratio `Psi4 / PySCF`

### C. Result comparison: C++ vs Python RESP core

Purpose:

- verify numerical parity for the migrated RESP algorithm

Use:

- identical ESP grid points
- identical ESP values
- identical equivalence constraints

This comparison must isolate the RESP solver from the QM backend so that any
delta is attributable to C++ migration rather than `Psi4` vs `PySCF`.

Compare:

- stage-1 charges
- final charges after stage-2
- max absolute delta
- RMS delta

Recommended first-cut tolerance:

- max absolute per-atom delta <= `1e-8`

If a different solver or ordering in C++ causes tiny floating-point drift,
document the exact reason before loosening tolerance.

### D. Speed comparison: C++ vs Python RESP core

Purpose:

- demonstrate the value of the C++ migration independent of QM runtime

Input:

- cached grid points
- cached ESP values
- cached topology-derived constraints

Measure:

- grid generation time
- matrix assembly time
- stage-1 solve
- stage-2 solve
- total RESP-core wall time

Success criterion:

- C++ core is measurably faster on the Mokda fixture
- no meaningful regression in charge agreement

## Suggested Benchmark Artifacts

Add a small benchmark script rather than relying only on pytest:

- `benchmarks/bench_resp_backends.py`

Suggested modes:

- `--backend pyscf`
- `--backend psi4`
- `--core python`
- `--core cpp`
- `--fixture mokda`
- `--fixture formamide`
- `--repeat N`

Suggested outputs:

- JSON summary for machine comparison
- human-readable console table

## Risk Register

### Risk 1. Psi4 ESP API details differ from PySCF

Impact:

- per-atom RESP charges may shift slightly even with the same basis

Mitigation:

- validate first on `only_esp=True`, `two_stage=False`
- compare raw ESP-on-grid values before comparing final RESP charges

### Risk 2. Geometry optimization changes coordinates differently

Impact:

- charge deltas may be dominated by geometry, not backend math

Mitigation:

- baseline migration comparisons should start with `opt=False`

### Risk 3. C++ migration introduces solver-order drift

Impact:

- tiny but noisy numerical differences

Mitigation:

- compare on cached ESP values
- keep deterministic matrix assembly order

### Risk 4. Benchmark noise hides real speed differences

Impact:

- misleading performance conclusions

Mitigation:

- separate QM backend timing from RESP-core timing
- use repeated runs and median reporting

## Initial Acceptance Criteria

The migration is complete when all of the following are true:

- `Assign.calculate_charge("RESP", ...)` remains stable
- `PySCF` remains the default backend
- `Psi4` backend passes the Mokda preview fixture
- C++ RESP core reproduces Python RESP core within defined tolerance
- benchmark evidence exists for:
  - `Psi4` vs `PySCF` result and speed comparison
  - C++ vs Python RESP-core result comparison
- docs and install guidance are updated, including Windows guidance for `Psi4`

## Immediate Next Steps

1. Add this plan file.
2. Refactor `resp.py` so the QM backend and RESP core are separable without changing behavior.
3. Add `Psi4` backend scaffolding behind an explicit non-default backend selection path.
4. Add fixture-driven comparison helpers for:
   - `PySCF` vs `Psi4`
   - Python core vs C++ core
5. Only after parity is demonstrated, move the RESP numerical core into C++.
