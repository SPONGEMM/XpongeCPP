# Xponge-origin `dev/assignment` Parity Migration Plan

This plan describes how to align `XpongeCPP` with the current
`Xponge-origin` `dev/assignment` branch.

The target is not a superficial file-by-file sync. The target is semantic
parity for assignment behavior, especially:

- bond-order refresh semantics
- AmberTools-aligned `GAFF` atom typing
- real `GAFF2` atom typing
- alternating-family orientation behavior
- sequence-sensitive RDKit-to-Assign behavior
- regression coverage for the new assignment edge cases

## Execution Rule

Before every implementation round for this migration, re-read this file first:

- [assignment-dev-branch-parity-plan.md](/mnt/data8t/Software/Xponge/Xponge-CPP/docs/migration/assignment-dev-branch-parity-plan.md)

This file is the source of truth for scope, sequencing, and acceptance.

## Scope

The migration source is the `Xponge-origin` branch:

- `dev/assignment`

The relevant upstream changes are concentrated in:

- `Xponge/assign/__init__.py`
- `Xponge/assign/bond_order.py`
- `Xponge/forcefield/amber/_alternating.py`
- `Xponge/forcefield/amber/gaff.py`
- `Xponge/forcefield/amber/gaff2.py`
- `Xponge/helper/rdkit.py`
- `Xponge/tools/unittests/test_1_assign.py`

## Goals

- Make `XpongeCPP` assignment behavior match `Xponge-origin/dev/assignment`
  for the changed `GAFF` and `GAFF2` cases.
- Keep the final implementation in the `XpongeCPP` C++ core where practical.
- Preserve existing `XpongeCPP` compatibility entry points such as
  `Assign.determine_atom_type(...)`.
- Add regression tests that compare against the new upstream semantics, not just
  against legacy `XpongeCPP` behavior.

## Non-goals

- Revert `XpongeCPP` back to a Python-only `GAFF` typing pipeline.
- Broaden this work into unrelated force-field, RESP, MCPB, or export topics.
- Rework the whole `Assign` API beyond what is needed for parity with
  `dev/assignment`.

## Current State in Xponge-origin

`Xponge-origin/dev/assignment` changed the assignment pipeline in a meaningful
way. The important parts are:

### 1. `Assign` now tracks bond insertion order

`Assign` has a `_bond_sequence` field. It is updated in:

- `add_bond`
- `delete_bond`
- `delete_atom`

This is not cosmetic. Later alternating-family post-processing depends on bond
sequence, not just the undirected bond graph.

### 2. Bond-order assignment explicitly invalidates derived state

After a successful bond-order update, upstream now does this:

- `assign.built = False`
- drop `_alternating_type_cache` if it exists

This forces ring markers, bond markers, aromatic markers, and alternating-family
adjustment state to be recomputed from the updated bond orders.

### 3. Alternating-family orientation is centralized

Upstream introduced:

- `Xponge/forcefield/amber/_alternating.py`

This file does a shared post-adjustment pass for:

- `cp/cq`
- `cc/cd`
- `ce/cf`
- `nc/nd`
- `ne/nf`
- `pc/pd`
- `pe/pf`

The upstream `GAFF` and `GAFF2` rules now rely on this post-pass instead of
encoding all direction decisions directly inside each rule.

### 4. `GAFF` was rewritten around the new post-pass

`Xponge/forcefield/amber/gaff.py` now:

- changes several primary rule predicates
- makes many secondary-family rules return `False` during the first pass
- uses `set_post_action(apply_amber_alternating_type_adjustment)` to orient the
  final atom types

### 5. `GAFF2` is now a real assignment implementation

`Xponge/forcefield/amber/gaff2.py` is no longer a data-only stub. It defines a
full typing rule set, including `GAFF2`-specific cases such as:

- `c5/c6`
- `ns/cs`
- special nitrogen behavior
- `s4` sulfoxide regression cases

It also uses the same alternating-family post-pass.

### 6. RDKit conversion behavior matters

`Xponge/helper/rdkit.py` upstream converts RDKit molecules to `Assign` while
preserving enough bond-order and traversal behavior to reproduce the new
sequence-sensitive regressions.

This matters for cases where:

- atom typing starts from SMILES or RDKit
- bond traversal order affects `ne/nf` versus `n2`

## Current State in XpongeCPP

`XpongeCPP` does not currently match the upstream architecture.

### 1. `GAFF` is implemented in C++

Current implementation:

- [cpp/assign/gaff.cpp](/mnt/data8t/Software/Xponge/Xponge-CPP/cpp/assign/gaff.cpp)

This file still encodes the older rule structure directly in C++ and does not
have the shared alternating-family post-pass used by upstream.

### 2. `GAFF2` is not really implemented yet

Current file:

- [src/XpongeCPP/forcefield/amber/gaff2.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/forcefield/amber/gaff2.py)

At the moment it only registers `gaff2.dat`. It is not a real `GAFF2` atom
typing pipeline.

### 3. `AssignRule` exists only as a compatibility path

Compatibility surface exists in:

- [src/XpongeCPP/assign/__init__.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/__init__.py)
- [src/XpongeCPP/_compat/assign.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/assign.py)

But `_compat/assign.py` explicitly keeps `gaff`, `gaff2`, and `sybyl` on the
core path instead of routing them through Python `AssignRule`.

### 4. Existing regression coverage is incomplete

`XpongeCPP` already has a small refresh regression:

- [tests/test_assign_bond_order_refresh.py](/mnt/data8t/Software/Xponge/Xponge-CPP/tests/test_assign_bond_order_refresh.py)

But it does not yet cover the full upstream `dev/assignment` regression set.

## Migration Principle

The migration target is:

- semantic parity with upstream
- implementation centered in the `XpongeCPP` C++ core

The migration target is not:

- copying upstream Python `GAFF` files verbatim into the compatibility layer

This means the correct strategy is:

1. port the upstream assignment semantics
2. keep `XpongeCPP` core ownership of `gaff` and `gaff2`
3. use Python only where `XpongeCPP` already depends on Python-facing glue,
   such as RDKit helpers or compatibility wrappers

## Target Architecture

### C++ core ownership

Core typing behavior should live in C++:

- `Assign` state tracking
- bond-order refresh invalidation
- `GAFF` typing
- `GAFF2` typing
- alternating-family post-adjustment

Suggested target implementation areas:

- `cpp/assign/assign.cpp`
- `cpp/core.hpp`
- `cpp/assign/bond_order.cpp`
- `cpp/assign/gaff.cpp`
- new `cpp/assign/gaff2.cpp`
- possibly a shared helper such as `cpp/assign/amber_alternating.cpp`

### Python compatibility ownership

Python should only own:

- legacy `AssignRule` compatibility
- RDKit conversion glue
- import-time force-field registration wrappers
- regression tests that exercise the public API

## Migration Phases

### Phase 0. Freeze baseline and comparison inputs

Deliverables:

- keep this plan in the repo
- identify the upstream regression cases to reproduce in `XpongeCPP`
- keep the current `XpongeCPP` assignment tests green before changes begin

Acceptance:

- plan file exists
- current assignment-related baseline tests still pass

### Phase 1. Port `Assign` state semantics

Goal:

Align `Assign` state transitions with upstream before touching `GAFF` logic.

Tasks:

- add `_bond_sequence` or a C++ equivalent to `Assign`
- update `add_bond` to record first-seen bond order insertion pairs
- update `delete_bond` to remove sequence entries
- update `delete_atom` to reindex and prune sequence entries
- align `built` invalidation semantics with upstream

Acceptance:

- state updates survive atom deletion and bond deletion
- benzene-style rebuild tests still pass
- the new sequence container is deterministic

### Phase 2. Port bond-order refresh invalidation

Goal:

Ensure successful bond-order assignment invalidates every derived state that must
be recomputed before typing.

Tasks:

- port the upstream refresh behavior from `bond_order.py`
- force marker/aromatic rebuild after bond-order changes
- invalidate any alternating-family cache if introduced in `XpongeCPP`

Acceptance:

- mol2 aromatic/ring-marker refresh regression passes
- repeated `determine_bond_order()` then `determine_atom_type()` is stable

### Phase 3. Implement shared alternating-family post-adjustment

Goal:

Create a shared AmberTools-aligned post-pass for alternating families.

Tasks:

- port the semantics of upstream `_alternating.py`
- implement:
  - primary-to-secondary normalization
  - bond-sequence-aware propagation
  - `cp/cq` special adjustment
  - sequence-sensitive nitroso demotion
  - specific `nc/nd` fallback behavior
- make the implementation reusable by both `gaff` and `gaff2`

Acceptance:

- dedicated `cp/cq` regression passes
- dedicated mixed `cc/cd` + `ce/cf` propagation regression passes
- nitroso `ne` and sequence-sensitive `n2` regressions pass

### Phase 4. Rewrite `GAFF` typing around the new post-pass

Goal:

Make `GAFF` behavior match upstream `dev/assignment`.

Tasks:

- update the core `GAFF` rule predicates to match upstream intent
- stop using the old direct local direction logic where upstream now uses a
  later post-pass
- call the alternating-family adjustment after the first typing pass

Acceptance:

- upstream `GAFF` regression cases pass
- `cc/cd` no longer degrade to `c2`
- `to_residuetype()` preserves corrected types

### Phase 5. Implement real `GAFF2` typing

Goal:

Replace the current `gaff2.dat` registration stub with a full `GAFF2` typing
implementation.

Tasks:

- add a real `GAFF2` rule set in the C++ core
- support `Assign.determine_atom_type("gaff2")`
- include `GAFF2`-specific rule families and reuse the alternating post-pass

Acceptance:

- benzene `gaff2` typing passes
- `c5/c6`, `ns/cs`, special nitrogen, and sulfoxide regressions pass

### Phase 6. Port RDKit helper semantics

Goal:

Ensure RDKit-derived assignments reproduce upstream sequence-sensitive behavior.

Tasks:

- add or update `XpongeCPP` RDKit helper glue so that bond insertion order is
  preserved into `Assign`
- match upstream aromatic/kekulization fallback semantics closely enough for the
  new regressions

Acceptance:

- RDKit/SMILES-origin regression cases reproduce upstream expected atom types

### Phase 7. Expand regression coverage

Goal:

Adopt the upstream regression intent, not just a subset.

Minimum test categories to port or adapt:

- bond-order refresh derived-state regression
- `GAFF` state and `to_residuetype()` regression
- cross-family alternating regression
- `GAFF2` special nitrogen cases
- `GAFF2` amide and thiocarbonyl cases
- `GAFF2` ring `c5/c6` cases
- `GAFF2` sulfoxide `s4` regression
- `cp/cq` pure aromatic orientation regression
- nitroso `ne` regression
- sequence-sensitive nitroso `n2` regression
- macrocycle kekule regression
- carbonyl ring `c` regression

Acceptance:

- `XpongeCPP` has dedicated tests for each upstream regression family

### Phase 8. Add origin-to-XpongeCPP parity checks

Goal:

Verify that `XpongeCPP` matches the current local upstream checkout for chosen
fixtures, not just hand-copied expectations.

Tasks:

- for a curated set of SMILES or mol2 fixtures, run both:
  - local `Xponge-origin`
  - local `XpongeCPP`
- compare final atom type sequences
- keep these checks isolated so optional dependency or environment differences
  do not destabilize unrelated CI jobs

Acceptance:

- parity fixtures show identical final atom-type sequences for the targeted
  `GAFF` and `GAFF2` cases

## Verification Strategy

Validation must happen at three layers.

### 1. Unit semantics

Verify local properties directly:

- `_bond_sequence` updates
- marker rebuild after bond-order changes
- alternating-family propagation behavior

### 2. Public API behavior

Verify through public entry points:

- `get_assignment_from_smiles`
- `get_assignment_from_mol2`
- `determine_atom_type("gaff")`
- `determine_atom_type("gaff2")`
- `to_residuetype()`

### 3. Cross-repo parity

Verify that local `XpongeCPP` matches local `Xponge-origin/dev/assignment`
results on a curated fixture set.

## Risks

### 1. C++ and Python ownership mismatch

Upstream now expresses `GAFF` and `GAFF2` in Python, while `XpongeCPP` owns
`GAFF` in C++ and does not yet own `GAFF2` there. This is the main architectural
translation risk.

### 2. Sequence-sensitive cases can fail silently

If `_bond_sequence` or its equivalent is missing or wrong, many ordinary tests
may still pass while the sequence-sensitive nitroso regressions drift.

### 3. Partial migration will produce misleading stability

If `GAFF` is updated without the shared post-pass, or if `GAFF2` is stubbed
instead of fully implemented, the result will look partially correct but still
diverge on exactly the cases that motivated `dev/assignment`.

## Completion Criteria

This migration is complete when all of the following are true:

- `XpongeCPP` has an explicit state model equivalent to upstream for
  bond-order-driven assignment refresh
- `GAFF` behavior matches `Xponge-origin/dev/assignment` on the migrated
  regression suite
- `GAFF2` is a real typing implementation, not only a parameter loader
- RDKit/SMILES-origin sequence-sensitive regressions are reproduced
- parity checks against local `Xponge-origin` pass for curated fixtures
- existing unrelated `XpongeCPP` workflows remain green

## First Implementation Order

The recommended execution order is:

1. `Assign` state semantics
2. bond-order refresh invalidation
3. shared alternating-family post-pass
4. `GAFF` rewrite
5. RDKit helper parity
6. real `GAFF2`
7. regression suite expansion
8. origin-to-XpongeCPP parity fixtures

This order minimizes false positives and avoids trying to debug `GAFF2` before
the lower-level state model is correct.
