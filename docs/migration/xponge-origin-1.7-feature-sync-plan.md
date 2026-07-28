# Xponge-origin 1.7 Feature Synchronization Plan

## Purpose and baseline

This is the execution plan for bringing `XpongeCPP` forward from the shared
explicit-periodic-box release to the relevant `Xponge-origin` 1.7 features.
It deliberately separates **integrating already-written XpongeCPP work** from
**implementing the remaining origin capabilities**.

Snapshot to re-check immediately before starting:

| Repository / ref | Commit | Meaning |
| --- | --- | --- |
| Xponge-origin `1.6b6` | `2ee3723` | Shared periodic-box baseline. |
| Xponge-origin `1.7b6` | `e60ef90` | Current source feature target. |
| XpongeCPP `master` | `eef6044` | Current release line (`v0.1.6`). |
| XpongeCPP `origin/codex/bundled-io` | `f28c159` | Existing native bundle-writer line. |
| Bundle branch merge base | `d14cb3f` | The integration point; it predates current `master`. |

The target is behavioral parity where it is meaningful in XpongeCPP's flat
C++ data model. It is not a source-level copy of Xponge-origin's nested Python
object model.

## Scope and ordering

```text
P0 baseline and integration branch
  -> P1 merge and validate native bundle writer
     -> P2 bundle read/conversion and numerical evidence
  -> P3 constrained RESP core and diagnostics
     -> P4 molecule-first metal-assignment contracts and transactional apply
        -> P5 manual_bonded overlays
           -> P6 QM/Hessian fitted-metal workflows and end-to-end releases
```

P2 and P3 may proceed in parallel after P1. P5 must not precede P4. P6 is
intentionally last because it combines the largest scientific and operational
validation surface.

## P0 — freeze evidence and create an integration branch

1. Refresh both remotes and record the five refs in the table above, including
   `git merge-base master origin/codex/bundled-io` and
   `git rev-list --left-right --count master...origin/codex/bundled-io`.
2. Start from a clean `master`; create a dedicated integration branch, for
   example `codex/integrate-bundled-io`.
3. Record the baseline outputs of the raw SPONGE saver and the current targeted
   test suite. Preserve generated artifacts outside the checkout.
4. Treat the existing MCPB draft and RESP migration plan as compatibility
   constraints; this plan supersedes neither their API guarantees nor their
   fixture strategy.

Exit criteria: the refs, divergence, test command versions, compiler and HDF5
versions are recorded in the PR/merge notes. No implementation starts from a
stale branch assumption.

## P1 — integrate `codex/bundled-io` as a merge, not a rewrite

The bundle branch is not fast-forwardable into current `master`: both sides
advanced after `d14cb3f`. Keep its commit history and merge it into the
integration branch with an explicit merge commit. Do not cherry-pick its
individual commits or reimplement `io_bundle.cpp`.

Resolve conflicts with these contracts:

- retain `master` behavior for raw text export and all existing public APIs;
- retain the branch's `format="bundle"`, `save_sponge_input_bundle`, native
  HDF5 writer, SHA256/UUID metadata, and atomic three-file publication;
- keep bundle an opt-in output (`raw` remains default);
- retain hard errors for unsupported bundle contents rather than emitting a
  partial or silently incompatible file;
- reconcile CMake, pixi, and Linux/macOS/Windows CI so HDF5 C and HighFive are
  discoverable on every wheel target.

The expected conflict hotspots are `CMakeLists.txt`, both wheel workflows,
`cpp/python/bindings_core.cpp`, `src/XpongeCPP/_compat/process.py`, public
exports, README/API documentation, and compatibility tests. New files such as
`cpp/io/io_bundle.cpp`, `cpp/io/sha256.hpp`, and
`tests/test_bundle_native.py` should remain the starting point unless a
conflict proves they are stale.

Validate before merging the integration branch back to `master`:

1. clean CMake configure/build with HDF5 and HighFive;
2. import smoke for raw and `format="bundle"` save paths;
3. `tests/test_bundle_native.py` and `tests/test_compat_surface.py`;
4. the branch parity helper against the sibling Xponge-origin checkout:
   `rtk uv run python scripts/validate_bundle_parity.py --xponge-root /mnt/data8t/Software/Xponge/Xponge-origin --work-dir <external-audit-dir>`;
5. the existing focused export/forcefield regressions affected by the merge;
6. Linux, macOS, and Windows wheel CI, including import smoke after packaging.

Exit criteria: bundle writer is on `master`, raw output remains unchanged for
the regression corpus, supported bundle schemas/hashes match the origin
fixtures, and all unsupported force modes fail explicitly.

## P2 — complete the bundle boundary without overstating parity

P1 provides native writing of topology, protocol sidecar, and H5MD restart;
it does not provide read-back or conversion parity. Add these in small,
independently reviewable changes:

1. typed topology/protocol/restart readers with schema-version validation;
2. a legacy-text to bundle converter and bundle-to-legacy converter for the
   documented supported-force subset;
3. an MDAnalysis adapter that reads these XpongeCPP bundle files, rather than
   relying on the unrelated optional analysis shim;
4. protocol objects beyond the empty/default sidecar only when their v2 schema
   contract and origin tests have been ported.

Required evidence:

- raw -> bundle -> parsed semantic equivalence;
- bundle -> legacy -> bundle equivalence for supported fixtures;
- corrupted schema, hash, lineage, atom-order, and path inputs fail clearly;
- at least one actual SPONGE numerical comparison (energy, forces or a short
  trajectory) before claiming runtime equivalence.

Do not silently support `fake_mass`, `fake_LJ`, `fake_charge`, arbitrary
listed forces, Stillinger-Weber, or EDIP merely by serializing them. Each needs
an explicit typed schema and an origin-parity test first.

## P3 — add constrained RESP as a standalone numerical capability

Source feature chain: origin `4fb8287` plus its charge-ledger, diagnostics,
preflight, worker-cleanup, permutation-stability, and open-shell hardening
commits through `1.7b5`.

Implement first in the existing RESP boundary:

- extend `src/XpongeCPP/assign/resp.py` and the C++ RESP numerical API to
  accept general linear constraints, equivalence groups, and target charge;
- use deterministic matrix assembly and report constraint residual, condition
  number, ESP-fit quality, and the applied constraint ledger;
- preserve current unconstrained RESP results and `Assign.calculate_charge`
  compatibility;
- keep QM backend execution separate from numerical fitting so PySCF/Psi4
  behavior can be compared with cached ESP inputs.

Acceptance tests must cover charge sum, equivalent atoms, arbitrary linear
constraints, infeasible/singular input, deterministic permutation behavior,
and Python-vs-C++ solver agreement on cached ESP data. This milestone can land
before metal assignment; metal-specific capped-model ledgers are a later
consumer of the same core.

## P4 — introduce molecule-first metal assignment

Origin's `bbe86ad`/`7efa65e` workflow is the design reference: a complete,
normally assigned parent `Molecule` receives a validated parameter overlay.
XpongeCPP should introduce a `metal_assignment` facade rather than force the
new model into the legacy `mcpb` implementation. `Xponge.MCPB()` remains a
compatibility entrypoint that adapts to the new request/result objects.

Build this in slices:

1. contracts and immutable artifacts: site mapping, electronic state,
   capped-model roles, topology/input hashes, and charge-accounting ledger;
2. selection/model construction using the flat `AtomId`/`ResidueId` model;
3. a parameterization plan that has no side effects on the parent molecule;
4. transactional `apply` that validates the plan, changes charges/bonded
   overlays/connectivity, rebuilds topology as required, and either publishes
   the whole result or leaves the parent unchanged;
5. saver integration and artifacts (`frcmod`, local model data, diagnostics).

The first supported production path requires explicit metal sites, atom
mapping, electronic state, and coordination edges. It must not infer a metal
coordination graph or silently cap a model. Validate Tier-1 metals from the
existing MCPB plan before expanding element coverage.

Acceptance: contract/artifact roundtrips; no mutation on validation failure;
atom/residue order preserved across successful apply; explicit metal links and
charges reach raw and bundle savers; and focused tests corresponding to
origin's contracts, base, charge, RESP, apply, and artifact suites.

## P5 — deliver `manual_bonded` as the first metal overlay mode

Port origin `2aa38d6` only after P4 provides stable topology/artifact identity.
Accept explicit reference bond lengths/angles and force constants, validate
term coverage and units, and tag every generated term as
`manual_bonded:explicit_reference_geometry`.

This mode must reject Hessian/RESP artifacts when they conflict with its
explicit-reference contract. It does not infer terms from distances and does
not claim QM fitting. Acceptance is a complete explicit geometry fixture,
coverage failures for omitted terms, deterministic overlay output, and raw/
bundle export of the applied parent molecule.

## P6 — fitted-metal workflow, release gate, and remaining parity

Only after P3-P5 are stable, add worker-isolated Hessian/Seminario fitting,
capped-model constrained RESP integration, robust open-shell SCF selection,
and lifecycle-safe worker cleanup. Treat every QM backend/basis/ECP combination
as a support matrix, not a generic promise.

Release only when all of the following hold:

- native bundle, constrained RESP, and each released metal mode have focused
  regression coverage and user-facing documentation;
- full non-QM test suite and platform wheel checks pass;
- affected QM tests pass or are explicitly gated with reproducible environment
  requirements;
- a real SPONGE numerical E2E result exists for supported raw and bundle paths;
- release notes state supported force types, metal modes, elements, QM
  backends, and known unsupported cases.

## Explicit non-goals for this program increment

- deleting the old MCPB API before its compatibility adapter is proven;
- automatic coordination/capping inference for arbitrary metalloproteins;
- full Amber MCPB.py, every special SPONGE force, or every origin CLI/analysis
  utility;
- claiming raw/bundle numerical equivalence from HDF5 schema or roundtrip tests
  alone;
- changing SPONGE runtime behavior itself.

