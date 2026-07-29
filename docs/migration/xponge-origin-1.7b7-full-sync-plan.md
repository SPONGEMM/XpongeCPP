# Xponge-origin 1.7b7 Full Synchronization Plan

## 1. Objective

Bring XpongeCPP from the last shared functional baseline to behavioral and
organizational parity with Xponge-origin `1.7b7`.

This is a full synchronization program, not three independent feature ports.
The target includes every functional change in Xponge-origin between:

| Project | Ref | Commit | Meaning |
| --- | --- | --- | --- |
| Xponge-origin | `1.6b6` | `2ee37236523f09fbadad9ac222e0a97c67de0760` | Last shared functional baseline: explicit periodic-box export |
| XpongeCPP | `v0.1.6` | `eef6044e5f908bd593429693cfcc4aae6c80c88c` | C++ implementation of the same periodic-box behavior |
| Xponge-origin | `1.7b7` | `2ef86d6d7be250c7b4ad57790eaa94d9b0a27880` | Synchronization target |

XpongeCPP `v0.1.5` is not the functional baseline. It was a corrective
dependency release corresponding to the earlier `1.6b5` feature state.

The final XpongeCPP design should resemble the origin module organization and
public contracts where practical. C++ implementations may replace Python
internals, but must not silently change public defaults, validation, metadata,
failure behavior, or artifact schemas.

## 2. Definition of done

The synchronization is complete only when all of the following are true:

1. Every non-release commit in
   `Xponge-origin/2ee3723..2ef86d6` is represented in a checked parity ledger
   as `implemented`, `not applicable`, or `intentional documented divergence`.
2. Bundle topology, protocol, restart, trajectory/output, conversion, and
   MDAnalysis behavior have origin-compatible typed contracts.
3. Generic RESP and QM orchestration expose the same constraints, defaults,
   progress events, SCF strategy/reference behavior, diagnostics, and metadata
   as origin. The C++ solver remains an explicit compatible extension.
4. `XpongeCPP.metal_assignment` owns the complete metal workflow: contracts,
   packages, structural/derived artifacts, base assignment, charge fitting,
   RESP/Hessian providers, force fitting, assigned-system materialization, and
   transactional molecule application.
5. No public or internal `MCPB` implementation remains:
   - no `XpongeCPP.MCPB` or `MCPB*` exports;
   - no `src/XpongeCPP/mcpb/`;
   - no `src/Xponge/mcpb/`;
   - no MCPB-only tests or documentation.
6. Origin's relevant test cohorts have been ported and pass against the C++
   data model.
7. Raw, native bundle, and bundle-to-legacy inputs produce matching numerical
   results in a real supported SPONGE runtime.
8. Linux, macOS, and Windows wheels build and pass their declared capability
   matrix.
9. Package version, CLI version, documentation, supported-feature table, and
   release notes agree.

## 3. Current branch assessment

The current `codex/integrate-bundled-io` branch is useful intermediate work,
not the final 1.7b7 port.

### Reuse after parity review

- native C++ bundle writer and SHA256/UUID publication;
- strict bundle reader and case discovery;
- legacy-to-bundle and bundle-to-legacy converters;
- basic MDAnalysis adapter;
- constrained RESP matrix preparation, diagnostics, permutation tests, and
  C++ KKT solver;
- C++ `Molecule` deep-copy replacement and molecule-local atom/LJ/bond/angle
  overrides;
- raw/bundle/SPONGE numerical test infrastructure.

### Replace or substantially restructure

- the simplified `metal_assignment` request/plan/result model;
- the standalone `prepare_manual_bonded_assignment` facade;
- the legacy MCPB-backed metal computation path;
- RESP orchestration that bypasses the origin-style QM scheduler;
- empty/default-only native protocol handling;
- reduced output/trajectory and MDAnalysis behavior.

Do not merge the branch to the release line until this plan's final gates pass.
Prefer forward corrective commits over rewriting already-visible history.

## 4. Workstream dependency graph

```text
P0 parity ledger and reproducible baseline
 ├─> P1 test isolation and build gates
 ├─> P2 bundle core contracts
 │    └─> P3 protocol, restart, output and analysis
 └─> P4 generic QM and RESP parity
      └─> P5 metal contracts and artifact packages
           └─> P6 base assignment and charge composition
                └─> P7 RESP/Hessian providers and worker lifecycle
                     └─> P8 force fitting and manual_bonded
                          └─> P9 assigned-system apply and MCPB removal
                               └─> P10 numerical, platform and release gates
```

P2/P3 and P4 may proceed in parallel after P0/P1. Metal work must not define
new private substitutes for contracts that already exist in origin.

## 5. Implementation phases

### P0 — freeze the parity ledger

Create `docs/migration/xponge-origin-1.7b7-parity-ledger.md` generated from:

```bash
git -C ../Xponge-origin log --reverse --no-merges \
  2ee37236523f09fbadad9ac222e0a97c67de0760..2ef86d6d7be250c7b4ad57790eaa94d9b0a27880
```

For each commit record:

- origin commit and release;
- feature/bug-fix behavior;
- origin source and tests;
- XpongeCPP owner files;
- status and intentional differences;
- verification command and evidence artifact.

Record compiler, HDF5/HighFive, Python, NumPy, h5py, MDAnalysis, QM backend,
CUDA, and SPONGE versions. Preserve baseline raw outputs outside the checkout.

Exit gate: no origin commit in the range is unclassified.

### P1 — make the test system trustworthy

The current process-global force-field registries make a monolithic pytest run
order-dependent. Fix this before interpreting full-suite results.

Implement one of:

1. complete C++ registry snapshot/restore exposed to a pytest fixture; or
2. fresh-interpreter grouping for every mutually incompatible force-field
   family.

The isolation must cover Amber templates and parameters, GAFF/GAFF2,
non-Amber registries, LJ combining rules, and Python module activation state.

Add:

- an order-randomized isolation test;
- ff14SB -> ff19SB and GAFF -> GAFF2 transition tests;
- a native C++ test target or an explicit Python ABI test for new C++ APIs;
- CI commands that run tests, rather than wheel import smoke alone.

Exit gate: two randomized full runs have the same pass/fail set, and no
failure is caused by state inherited from an earlier test.

### P2 — synchronize bundle core contracts and native writing

Port the origin `io_bundle` organization while retaining the native C++ writer:

- schema/constants/contracts;
- canonical topology, force-field, protocol, restart-state, content, and
  lineage hashing;
- topology/state/trajectory parser and exporter contracts;
- molecule adapter boundary;
- typed HDF5 writer utilities;
- atomic three-file publication and rollback;
- CLI/case discovery and manifest contracts.

The C++ writer should produce the same logical contract directly. Do not
generate raw text and reparse it as the implementation.

Port origin bundle contract and saver tests, then compare:

- schema/version fields;
- UUID and parent/source lineage;
- atom order and topology hashes;
- restart finalized/load-policy behavior;
- corrupt/truncated publication rejection;
- unsupported force errors.

Exit gate: origin and C++ writers pass the same canonical fixture validator for
the supported force matrix.

### P3 — synchronize protocol, restart, output, conversion, and analysis

Implement the full typed protocol surface used by origin `1.7b7`:

- collective variables;
- constraints and restraints;
- metadynamics;
- steering;
- SITS;
- hard/soft walls;
- protocol serialization and hashing.

Then synchronize:

- restart state hash and completion metadata;
- native output/H5MD writer and legacy output ingestion;
- topology/state/trajectory parsers and exporters;
- multi-frame and walker/stream semantics;
- units and layout handling in MDAnalysis;
- both conversion directions with an explicit capability matrix.

Do not claim support for a force merely because it can be stored as an opaque
array. Unsupported typed contents must fail with a stable error.

Exit gates:

- every protocol dataclass round-trips deterministically;
- raw -> bundle -> legacy -> bundle semantic parity;
- topology plus multi-frame MDAnalysis parity;
- corrupted schema/hash/lineage/unit/layout fixtures fail as origin does.

### P4 — synchronize generic QM and RESP

Keep the existing constrained solver work, but restore the full origin public
behavior around it.

#### Public RESP parity

Match origin `resp_fit` parameters and defaults:

- `progress_callback`;
- `scf_strategy`;
- `scf_reference`;
- ESP memory/chunk/safety options;
- constraint matrix/targets;
- metadata and diagnostics return shapes.

`core=None` must retain origin-compatible behavior. Keep `core="cpp"` as an
explicit XpongeCPP extension until a deliberate compatibility-version decision
changes the default.

#### Scheduler and backend parity

Synchronize:

- `QMRunOptions` memory limit, SCF tolerance, maximum cycles, strategy, and
  reference;
- `resolve_scf_reference(auto/rhf/rohf/uhf)`;
- `SCFResult.reference`;
- Hessian convergence, energy, and reference metadata;
- scheduler-based SCF/ESP/optimization/Hessian execution;
- PySCF RHF/ROHF/UHF and direct/density-fit/newton selection;
- Psi4 reference handling;
- exact physical conversion constants.

Generic RESP must:

- reject unconverged SCF;
- emit ordered progress phases;
- report timings, backend, requested/resolved reference, energy, grid count,
  and ESP diagnostics;
- preserve all constraints through initial, stage-1, and stage-2 solves;
- preflight invalid constraints before starting QM.

Port origin `test_1_charge.py` and `test_1_resp_linear_constraints.py`,
including open-shell, callback, metadata, second-stage hydrogen, constraint,
and permutation cases. Add malformed-input tests for the direct C++ solver.

Exit gate: origin-compatible Python core tests pass; `core="cpp"` agrees on
cached ESP fixtures within the declared tolerance.

### P5 — replace the simplified metal contracts with origin contracts

Mirror the origin package organization:

- `contracts.py`;
- `artifacts.py`;
- `input.py`;
- `apply_input.py` where applicable;
- deterministic JSON serialization and schema validation.

Implement:

- prepared atoms, residues, components, bonds, and links;
- chemical topology and residue-partition proofs;
- electronic-state and charge contracts;
- provider capability snapshots and projections;
- structural artifacts;
- derived small/large models and closed atom mappings;
- `PreparedArtifactPackage`;
- `MetalAssignmentPackage`;
- `ParameterizationResult`.

Adapt storage to flat `AtomId`/`ResidueId`, but preserve external IDs, hashes,
units, graph revisions, mappings, and validation errors.

Include the `1.7b7` solvent-ion versus metal-center classification correction;
do not classify every single-atom ion as a fitted metal center.

Exit gate: all contract/package/artifact serializer and stale-hash tests ported
from origin pass byte-deterministically.

### P6 — base assignment, provider projection, and charge composition

Implement origin-equivalent modules:

- base and standard force-field assignment;
- template and cross-provider assignment;
- base composition;
- nonbonded metal-ion assignment;
- metal overlay construction;
- base charge inputs and providers;
- partial-charge artifact composition and projection;
- assigned-system artifact construction.

Worker/provider outputs must be immutable and hash-closed. Provider
capabilities, water model, basis/ECP, force protocol, source package hash, and
application audit must remain visible.

Exit gate: origin base/base-charge/nonbonded tests pass, including missing
coverage, provider conflict, solvent-ion classification, and deterministic
projection failures.

### P7 — RESP/Hessian providers and worker lifecycle

Build the metal-local provider layer on P4 rather than reusing MCPB helpers:

- `RespFitInput`, constraints, equivalence groups, and serialization;
- constrained local-model RESP worker;
- `ModelChargeArtifact` and projection to the parent;
- `HessianFitInput` and Hessian worker;
- backend/basis/ECP/reference capability validation;
- worker protocol versioning;
- timeout, cancellation, process-group termination, stderr preservation, and
  guaranteed reaping.

The charge ledger must explicitly account for parent atoms, retained model
atoms, caps, boundary redistribution, and final total charge.

Exit gate: origin metal RESP/Hessian test cohorts pass with cached providers;
real backend tests run only on declared platform/backend combinations.

### P8 — synchronize force fitting and `manual_bonded`

Implement under the origin-style modules:

- `BondedFitInput`;
- atom and bonded metal parameter specifications;
- Seminario terms;
- empirical registry and applicability metadata;
- `manual_bonded_terms`;
- bonded-fit composition.

Retain useful molecule-local C++ bond/angle/LJ override machinery, but make it
an application backend for origin-compatible `ParameterizationResult`, not a
parallel public contract.

`manual_bonded` must:

- consume an immutable explicit-reference artifact;
- require frozen current geometry and explicit units;
- cover every active coordination bond and required donor-pair angle exactly;
- reject RESP/Hessian artifacts;
- use `manual_bonded:explicit_reference_geometry` provenance;
- never infer equilibrium values from an unconfirmed current structure.

Exit gate: manual, empirical, and Seminario fixtures produce deterministic
terms and raw/bundle exports; missing/duplicate/extra coverage fails clearly.

### P9 — assigned-system application and complete MCPB removal

Implement origin-compatible:

- assigned-system artifacts and force realization protocol;
- isolated materialization worker where required;
- molecule-first `assign`, `parameterize`, and `apply`;
- topology fingerprint preservation;
- registry snapshot/rollback;
- copy-first application and one-step `inplace=True` commit.

Then remove, in the same milestone:

- `src/XpongeCPP/mcpb/`;
- `src/Xponge/mcpb/`;
- all `MCPB`, `MCPBIonInfo`, `MCPBRequest`, `MCPBResult`,
  `MCPBSelection`, and `MCPBLocalModel` exports;
- MCPB-only documentation and tests.

Do not retain a deprecated wrapper unless the user explicitly reverses this
decision. Rewrite useful MCPB fixtures against `metal_assignment`; do not
delete scientific coverage with the old namespace.

Exit gate:

```bash
rg -n 'MCPB|XpongeCPP\\.mcpb|Xponge\\.mcpb' \
  src tests docs README.md
```

returns no product/API references, apart from an intentional migration note.

### P10 — full E2E, platform, and release gate

Run mandatory, non-skipped numerical gates on a dedicated SPONGE runner:

1. ordinary force field: raw versus native bundle/converted legacy;
2. constrained RESP fixture;
3. nonbonded metal-ion fixture;
4. `manual_bonded` fixture;
5. one fitted Seminario/Hessian fixture on a declared backend.

Compare all zeroth-frame mdout columns, then at least one short deterministic
trajectory or minimization for energies, forces, and coordinates.

CI matrix:

- Linux x86_64 and aarch64;
- macOS Intel and arm64;
- Windows x86_64;
- every supported CPython version;
- HDF5/HighFive discovery;
- import, focused pytest, wheel install, and artifact smoke;
- platform-specific QM capability expectations.

Also add:

- `python -m XpongeCPP -v` equals `XpongeCPP.__version__`;
- plain import has no stdout/stderr;
- optional SPONGE tools use environment variables or `PATH`, never developer
  absolute paths;
- `twine check`, wheel metadata/dependency checks, and reproducible hashes.

Because MCPB removal is a breaking API change, use XpongeCPP `v0.2.0` unless a
different versioning decision is made explicitly. Record
`Xponge-origin compatibility target: 1.7b7` in package metadata and release
notes rather than reusing origin's tag scheme.

Exit gate: all required CI and numerical jobs pass; no required job is skipped.

## 6. Recommended review/commit boundaries

Use small reviewable commits or PRs in this order:

1. parity ledger and deterministic test isolation;
2. bundle contracts/hash/writer;
3. protocol/restart/output/analysis;
4. generic QM models/scheduler/backends;
5. generic RESP API/tests and C++ parity;
6. metal contracts/artifact packages;
7. base assignment and charge composition;
8. RESP/Hessian providers and worker runtime;
9. bonded fit/manual/empirical/Seminario;
10. assigned-system apply and C++ overlay backend;
11. MCPB deletion and public API cleanup;
12. numerical/platform/release hardening.

Each commit must include its focused tests and update the parity ledger. Avoid
a final bulk test-port commit disconnected from the implementation it verifies.

## 7. Scope controls

The following do not count as parity:

- retaining MCPB and wrapping only its final mutation in a transaction;
- matching class names while omitting serialization, hashes, provenance, or
  failure behavior;
- claiming bundle support from empty protocol sidecars;
- claiming numerical parity from schema or round-trip tests alone;
- silently switching RESP defaults because a C++ implementation exists;
- accepting a backend/basis/ECP combination without a tested support-matrix
  entry;
- declaring the full suite green by excluding order-dependent failures.

Any intentional divergence must be recorded in the parity ledger with a user-
visible reason, compatibility effect, and dedicated test.
