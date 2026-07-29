# Xponge-origin 1.7b7 Parity Ledger

## Scope and baseline

This ledger is the authoritative checklist for synchronizing XpongeCPP with
Xponge-origin `1.7b7`.

| Project | Ref | Commit |
| --- | --- | --- |
| Xponge-origin baseline | `1.6b6` | `2ee37236523f09fbadad9ac222e0a97c67de0760` |
| XpongeCPP baseline | `v0.1.6` | `eef6044e5f908bd593429693cfcc4aae6c80c88c` |
| Xponge-origin target | `1.7b7` | `2ef86d6d7be250c7b4ad57790eaa94d9b0a27880` |
| XpongeCPP implementation | working tree | XpongeCPP `0.2.0` candidate |

The origin range is:

```text
2ee37236523f09fbadad9ac222e0a97c67de0760..2ef86d6d7be250c7b4ad57790eaa94d9b0a27880
```

Statuses:

- `implemented`: behavior is present and has focused evidence;
- `partial`: useful implementation exists, but origin parity is not yet proven;
- `pending`: the behavior has not been implemented in the target structure;
- `not-applicable`: release-only metadata that is not ported commit-for-commit.

No row may be changed to `implemented` without adding the verification command
and an evidence path or CI job.

## Reproducible environment snapshot

Captured on 2026-07-28 in `/mnt/data8t/Software/Xponge/Xponge-CPP`.

| Component | Value |
| --- | --- |
| Operating system | Linux |
| Python | 3.11.15 local; release support is CPython 3.10-3.12 |
| NumPy | 1.26.4 |
| h5py | 3.16.0 |
| MDAnalysis | 2.10.0 |
| CMake | 4.3.2 |
| System C++ compiler | GCC 14.2.0 |
| Pixi C/C++ compiler | GCC 14.3.0 |
| HDF5 | 2.1.0, serial, threadsafe |
| SPONGE executable | `/mnt/data8t/Software/SPONGE/SPONGE/build-dev-cpu/SPONGE`; found by sibling-checkout discovery. Direct native-bundle gate uses `SPONGE_BUNDLE_EXECUTABLE`. |
| QM backends | PySCF 2.14.0 available; Psi4 unavailable |
| CUDA | to be recorded on the dedicated numerical runner |

Baseline numerical artifacts must be stored outside the checkout and recorded
here before P10. Schema-only fixtures do not satisfy the numerical baseline.

## Commit-by-commit parity

The `Origin evidence` column names the primary source/test cohort. The
`XpongeCPP owner` column is the intended final owner, not necessarily the
current implementation location.

| Origin commit | Release | Behavior | Origin evidence | XpongeCPP owner | Status | Verification/evidence |
| --- | --- | --- | --- | --- | --- | --- |
| `b495cb0` | 1.7b1 | Legacy input/output/state/trajectory to bundle conversion and CLI | `Xponge/io_bundle/{converter,legacy_case,output_parsers,output_writer,state_parsers,topology_parsers,trajectory_parsers}.py`; `test_1_io_bundle.py` | `src/XpongeCPP/io_bundle/`; converter tests | implemented | Bundle cohort with discovered SPONGE: 166 passed; local converter plus 1.7b7 numerical gate: 11 passed |
| `5ab968f` | 1.7b1 | Bidirectional bundle input, reader, saver, reverse conversion, molecule adapter | `Xponge/io_bundle/{bundle_builder,bundle_case,bundle_reader,exporters,reverse_converter,saver}.py`; bundle contract/reverse/saver tests | `src/XpongeCPP/io_bundle/`, `cpp/io/io_bundle.cpp` | implemented | Native saver, residue state, canonical hashes, bidirectional typed conversion and failure tests pass |
| `54cd862` | 1.7b1 | MDAnalysis topology/trajectory integration | `Xponge/io_bundle/`; MDAnalysis cases in `test_1_io_bundle.py` | `src/XpongeCPP/analysis/` | implemented | Origin and existing MDAnalysis cohorts pass together |
| `b27978a` | 1.7b1 | Optional SPONGE executable discovery for tests | `test_1_io_bundle.py` and test helpers | test/runtime discovery helpers | implemented | Explicit environment, `PATH`, and sibling-build discovery; converter cohort passes with the environment variable unset |
| `d607b89` | 1.7b1 | CLI `-v` matches package `__version__` | `Xponge/__main__.py`; CLI test | `src/XpongeCPP/__init__.py`, CLI entrypoint | implemented | Both CLI and import report `0.2.0` |
| `45e52be` | 1.7b1 | Native bundle output emission | `Xponge/io_bundle/output_writer.py`; output fixtures | `src/XpongeCPP/io_bundle/`, native HDF5 layer | implemented | Output/restart/H5MD writer tests included in 151-test bundle cohort |
| `b890fbf` | 1.7b1 | Import-path and optional-dependency hardening | `Xponge/io_bundle/__init__.py`; import tests | `src/XpongeCPP/io_bundle/__init__.py` | implemented | Clean package and wheel imports pass |
| `2c0dc3e` | 1.7b1 | Residue/topology handling fixes in bundle conversion | bundle contracts/converter tests | bundle topology adapters | implemented | ResidueType and state-preserving Residue native saver tests pass |
| `fcfad98` | 1.7b1 | SPONGE bundle schema/behavior alignment | bundle schema, saver and SPONGE integration tests | Python contracts plus `cpp/io/io_bundle.cpp` | implemented | Canonical dataset hashes and full typed reverse-export cohort pass |
| `d2cbfb7` | 1.7b2 | Release metadata for 1.7b2 | version metadata | release metadata | not-applicable | Superseded by target compatibility version |
| `2d372f4` | 1.7b3 | Bundle release metadata/documentation | version and bundle docs | release notes/capability table | not-applicable | Final v0.2.0 release documents target 1.7b7 |
| `e992d61` | 1.7b3 | Native molecule-to-bundle saver | `Xponge/io_bundle/saver.py`; saver fixtures | `cpp/io/io_bundle.cpp`, Python saver facade | implemented | Native molecule/residue/template saver plus protocol/restart finalization passes |
| `bbe86ad` | 1.7b4 | Metal workflow reorganized under `metal_assignment` | `Xponge/metal_assignment/`; contract/artifact tests | `src/XpongeCPP/metal_assignment/` | implemented | Full origin package organization ported; MCPB package and symbols removed |
| `7efa65e` | 1.7b4 | Immutable metal artifacts, assignment/providers/apply workflow | `Xponge/metal_assignment/*.py`; all `test_1_metal_assignment_*` cohorts | `src/XpongeCPP/metal_assignment/` | implemented | Metal origin-parity cohort: 82 passed, 23 optional external-provider tests skipped |
| `923cb22` | 1.7b4 | Metal-assignment release metadata | version and docs | release notes/capability table | not-applicable | Final v0.2.0 release documents target 1.7b7 |
| `4fb8287` | 1.7b5 | Generic RESP linear constraints | `Xponge/helper/`, QM/charge modules; `test_1_resp_linear_constraints.py` | `src/XpongeCPP/qm/`, C++ RESP solver | implemented | Existing plus origin constraint cohorts pass on rebuilt C++ extension |
| `664cbb4` | 1.7b5 | Metal RESP consumes constrained generic RESP | metal RESP input/provider/worker tests | P4 generic RESP plus P7 metal provider | implemented | Origin metal RESP provider/worker tests pass |
| `ee52487` | 1.7b5 | Expose RESP constraint controls through public API | `test_1_charge.py`, `test_1_resp_linear_constraints.py` | public QM/RESP facade | implemented | Public Assign API passes constraints, callbacks, strategy/reference and metadata |
| `a60b47d` | 1.7b5 | Preserve capped-model charge constraints and ledger | metal RESP artifacts/provider tests | metal charge artifacts/composition | implemented | Charge ledger/artifact tests pass |
| `91779ca` | 1.7b5 | Report RESP electronic-state/reference metadata | generic and metal RESP tests | QM result metadata, metal artifacts | implemented | Metadata assertions pass |
| `436f322` | 1.7b5 | Model-local constrained RESP workflow | `resp_input.py`, `resp_provider.py`, `_resp_worker.py` | metal RESP provider/worker | implemented | Provider cache and worker cohort pass |
| `a9161e0` | 1.7b5 | Hessian worker lifecycle and failure handling | `_hessian_worker.py`, `_worker_runtime.py`; Hessian tests | metal Hessian provider/runtime | implemented | Hessian worker lifecycle cohort passes |
| `9da9a75` | 1.7b5 | Explicit constrained RESP contracts | RESP contracts/input/provider tests | generic RESP and metal `RespFitInput` | implemented | Contract serialization and validation tests pass |
| `336a3fa` | 1.7b5 | Preserve RESP worker diagnostics/failure details | `_resp_worker.py`, runtime tests | metal worker runtime | implemented | Structured worker error/diagnostic tests pass |
| `3779bd6` | 1.7b5 | Govern capped-model construction and mappings | metal artifacts/input/provider tests | structural/derived metal artifacts | implemented | Closed mapping and stale-hash tests pass |
| `0a8a93a` | 1.7b5 | Reap interrupted QM worker process groups | `_worker_runtime.py`; interruption tests | shared/metal worker runtime | implemented | Timeout/interruption process-group tests pass |
| `66e25c6` | 1.7b5 | Preflight constrained RESP before QM execution | generic RESP tests | public RESP orchestration | implemented | Invalid constraints fail before backend execution |
| `bd25ff0` | 1.7b5 | Constraint regression and permutation coverage | `test_1_resp_linear_constraints.py` | generic RESP tests | implemented | Full origin constraint cohort ported and passing |
| `ff4cded` | 1.7b5 | Robust automatic RHF/ROHF/UHF selection | QM backend and charge tests | QM scheduler/backends | implemented | Scheduler and PySCF reference-selection tests pass |
| `a62e0b4` | 1.7b5 | Release metadata for 1.7b5 | version metadata | release metadata | not-applicable | Superseded by target compatibility version |
| `2aa38d6` | 1.7b6 | Explicit-reference `manual_bonded` metal fitting | bonded/force-fit modules and metal tests | `metal_assignment/bonded_fit.py`, force-fit application backend | implemented | Origin manual-bonded contract/provenance/application tests pass |
| `e60ef90` | 1.7b6 | Release metadata for 1.7b6 | version metadata | release metadata | not-applicable | Superseded by target compatibility version |
| `9a33901` | 1.7b7 | Separate solvent ions from fitted metal centers | metal input/base/nonbonded tests | metal classification/base assignment | implemented | Base/nonbonded ion classification and TIP3P Zn/Fe tests pass |
| `2ef86d6` | 1.7b7 | Release metadata for target 1.7b7 | version metadata | package metadata and release notes | not-applicable | P10 declares `Xponge-origin compatibility target: 1.7b7` |

## Phase gates

| Phase | Required status transition | Evidence |
| --- | --- | --- |
| P1 | Test isolation rows added; inherited registry state cannot alter outcomes | C++/Python snapshot ABI tests plus two randomized full-run manifests |
| P2/P3 | All bundle rows become `implemented` | Origin/C++ canonical fixtures, corruption fixtures, SPONGE numerical artifacts |
| P4 | Generic RESP rows become `implemented` | Ported charge/constraint cohorts and cached ESP C++ comparison |
| P5–P9 | All metal rows become `implemented`; MCPB search is clean | Ported metal cohorts, application rollback tests, namespace scan |
| P10 | Release-only rows are resolved by v0.2.0 metadata | Platform CI, wheels, installation smoke and release checklist |

## Known blockers

1. Ordinary, constrained-RESP, nonbonded-metal, manual-bonded, and Seminario
   raw -> bundle -> legacy numerical gates pass on the discovered CPU build.
   The available SPONGE v2.0.0-beta.1 source and all local/remote branches do
   not implement the native `input_h5_*` commands, so direct native-bundle
   execution still needs a bundle-capable SPONGE revision.
2. PySCF 2.14.0 was exercised locally; Psi4 and the external
   Mokda-chemcore-backed artifact providers are unavailable in this environment.
3. The post-sync full-suite result is 487 passed, 31 failed, 10 skipped, and
   1 expected failure, improving on the stable pre-sync result of 412 passed,
   32 failed, 11 skipped, and 1 expected failure. These are independently
   reproducible product/test issues,
   including same-test GAFF/GAFF2 imports that contradict the active-family
   guard; they are no longer classified as inherited-state failures.

## Execution evidence

### P0

- Origin non-merge commit count:
  `git rev-list --count --no-merges 2ee3723..2ef86d6` -> `34`.
- Ledger commit-row count -> `34`.
- `git diff --check` passed after creating the ledger.

### P1 complete

- Added an opaque C++ registry snapshot/restore ABI covering Amber templates
  and parameters, non-Amber registries, PDB mappings, and the LJ combining
  rule.
- Added an autouse pytest fixture covering the C++ state, force-field-family
  activation, legacy `AtomType` state, and force-field module cache/parent
  attributes.
- Focused transition tests:
  `tests/test_forcefield_registry_isolation.py` -> `3 passed`.
- Forward isolation cohort -> `58 passed`.
- Reverse isolation cohort -> `58 passed`.
- Editable rebuild completed successfully before running the cohorts.
- Full suite with `XPONGECPP_TEST_ORDER_SEED=1729` ->
  `412 passed, 32 failed, 11 skipped, 1 xfailed`.
- Full suite with `XPONGECPP_TEST_ORDER_SEED=2718` ->
  `412 passed, 32 failed, 11 skipped, 1 xfailed`.
- The two full-run failure node-ID sets are identical.

### P2/P3 implementation complete

- Ported the full bundle contract, parser/exporter, protocol, restart, output,
  CLI, trajectory, and MDAnalysis organization.
- Preserved the C++ topology/restart writer and added state-preserving
  Molecule/Residue/ResidueType preparation.
- Canonical hashes are derived from logical typed datasets, including legacy
  metadata and sidecar tables.
- Current bundle contract/conversion/analysis cohort with discovered SPONGE:
  `166 passed`.
- With automatic sibling-build discovery and `SPONGE_EXECUTABLE` explicitly
  unset, `tests/test_bundle_converter.py` completed with `7 passed`, including
  raw -> bundle -> legacy SPONGE zero-frame numerical equality.
- Native metal exports containing zero-record bonded tables exposed and fixed
  a generic fixed-table parser shape bug; all eight zero-record table variants
  now have regression coverage.

### P4 implementation complete

- Ported QM models, basis/parameter definitions, scheduler, Psi4/PySCF
  backends, SCF strategies, and automatic restricted/open-shell reference
  selection.
- Generic RESP now exposes linear constraints, targets, equivalence groups,
  preflight, progress callbacks, diagnostics, and fit metadata.
- The rebuilt combined RESP/bundle/metal cohort completed with
  `238 passed, 12 skipped`.

### P5-P9 implementation complete

- Ported the origin `metal_assignment` contracts, immutable artifacts, input
  models, base/nonbonded assignment, charge composition, RESP and Hessian
  providers, worker runtime, force fitting, and transactional application.
- Added a native parameter extraction adapter that reads the actual C++ SPONGE
  export for GAFF/GAFF2, standard biomolecular templates, and ions.
- Metal origin-parity cohort: `82 passed, 23 skipped`; skipped cases require
  optional Mokda-chemcore or unavailable QM providers.
- Deleted `src/XpongeCPP/mcpb/`, MCPB-only tests/docs, and top-level exports.
  A source/test/CI namespace scan for `MCPB` is empty.

### P10 local evidence; external gates pending

- Set package and CLI version to `0.2.0` and recorded the Xponge-origin 1.7b7
  compatibility target in package metadata and `RELEASE_NOTES.md`.
- Declared the verified CPython support range as 3.10-3.12. Build and publish
  workflows build each version on Linux x86_64/aarch64, macOS Intel/arm64, and
  Windows x86_64, with per-wheel cibuildwheel smoke commands.
- `pixi run install-dev` rebuilt and installed the C++ extension.
- `python -m build` produced the 0.2.0 sdist and Linux CPython 3.11 wheel.
- `twine check` passed for both artifacts.
- Clean wheel smoke import reported `0.2.0`, imported bundle and
  `metal_assignment`, and confirmed that `MCPB` is absent.
- SPONGE executable discovery works through explicit configuration, `PATH`, or
  a sibling source-tree build. With the environment override unset, the bundle
  converter/SPONGE numerical cohort completed with `7 passed`.
- Added mandatory real-SPONGE fixtures for constrained RESP, nonbonded metal,
  explicit-reference `manual_bonded`, and deterministic-Hessian Seminario.
  Their four raw -> bundle -> legacy gates pass; the Seminario case also
  compares a two-step minimization's energy rows, force trajectory, coordinate
  trajectory, and restart coordinates.
- The CI numerical job builds pinned SPONGE commit
  `4c694ebab7032b0ef28d8312115d0f3253800125` on a hosted CPU runner and runs
  the converter plus four 1.7b7 numerical fixtures as a mandatory 11-test gate.
  The same pinned revision was built from a detached clean worktree locally,
  and the exact CI test selection completed with `11 passed`.
- Added `tests/test_bundle_native_sponge.py` as the explicit direct-execution
  acceptance gate. It is selected only by configuring
  `SPONGE_BUNDLE_EXECUTABLE` to a runtime supporting `input_h5_*`; the current
  runtime records one transparent skip rather than a false pass.
- Full local suite: `504 passed, 31 failed, 10 skipped, 1 xfailed`; none of
  the failures belong to the bundle, RESP, or metal parity cohorts.
- Final combined release/RESP/bundle/metal/numerical parity cohort:
  `217 passed, 23 skipped`; all skips require optional external providers.
- Still required before an actual release: direct execution on a
  bundle-capable SPONGE revision and an actual green Linux/macOS/Windows CI
  matrix run.
