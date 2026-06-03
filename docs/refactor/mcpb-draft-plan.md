# Xponge MCPB Draft Plan

This document is an implementation draft for an `Xponge.MCPB()` workflow in
`XpongeCPP`. The goal is to provide a metal-center parameterization path that
is conceptually similar to AmberTools `MCPB.py`, while fitting the current
`XpongeCPP` architecture:

- flat `Molecule` / `Residue` / `Atom` C++ storage
- `Assign` for local chemistry and charge workflows
- shared `qm/` subsystem for SCF / ESP / optimization
- C++ RESP numerical core

This draft is intentionally pragmatic. It focuses on a first implementation
path that can grow into fuller MCPB-like support later.

## Goals

- Add an `Xponge.MCPB()` user-facing workflow for metal centers.
- Accept a loaded `Molecule` as the primary input rather than requiring a
  whole-system `Assign`.
- Support explicit user selection of metal ion atoms and nearby coordination
  environment.
- Accept all periodic-table elements at the API/input layer, while separating
  validated support tiers from broad input acceptance.
- Make the post-MCPB parent `Molecule` directly exportable through
  `Save_SPONGE_Input(...)` / `save_sponge_input(...)`.
- Generate and integrate:
  - metal-centered bonded parameters
  - locally refit ligand/residue charges
  - explicit connectivity information for metal bonds
  - metal ion typing / template support needed by topology and SPONGE export
- Preserve auxiliary artifacts such as `frcmod`, local model files, and
  metadata for inspection and debugging.
- Reuse the shared QM layer instead of building a second RESP/QM stack.

## Non-goals For The First Version

- Full parity with every `MCPB.py` step and switch on day one.
- Fully automatic coordination inference for arbitrary metalloproteins.
- All bonded and nonbonded MCPB variants in the first release.
- Uniform day-one validation, parameter availability, and benchmark coverage
  for every element in the periodic table.

## Why `Molecule` Should Be The Main Input

For a whole metalloprotein or metal-containing complex, `Molecule` is the safer
starting point than `Assign`:

- `load_pdb` / `load_mol2` can usually read metal-containing systems into a
  `Molecule` without forcing bond-order perception.
- `Assign` currently performs connectivity and bond-order logic that is mainly
  organic-small-molecule oriented.
- Many transition-metal centers can fail in `Assign.determine_bond_order(...)`
  with `No bond-order atom type for atom #...`.

Therefore the MCPB workflow should:

1. read the full system as a `Molecule`
2. let the user explicitly choose the metal ion(s)
3. build small and large local models from the selected region
4. use `Assign` only for local fragments when appropriate

## Proposed User-facing API

The long-term public entrypoint should look like:

```python
result = Xponge.MCPB(
    molecule,
    ion_ids=[1523],
    ion_info=[
        {
            "atom_id": 1523,
            "element": "Zn",
            "formal_charge": 2,
            "spin": 1,
            "resname": "ZN",
            "atom_name": "ZN",
        }
    ],
    method="seminario",
    model="bonded",
    cutoff=2.8,
    bonded_pairs=[(1523, 1472), (1523, 1489), (1523, 1501), (1523, 1510)],
    additional_residue_ids=[84, 109],
    charge_mode="resp",
    qm_backend="pyscf",
    basis="6-31g*",
    scale_factor=1.0,
)

result.molecule.Save_SPONGE_Input("case", dirname="out")
```

## Proposed Minimal API Contract

### Required inputs

- `molecule`
  - a loaded `XpongeCPP.Molecule`
- `ion_ids`
  - atom ids of the target metal centers

### Strongly recommended explicit inputs

- `ion_info`
  - list of dicts containing at least:
    - `atom_id`
    - `element`
    - `formal_charge`
    - `spin`
- `bonded_pairs`
  - explicit metal-ligand bonded pairs

### Optional inputs

- `method`
  - force-constant generation method
- `model`
  - `bonded` or `nonbonded`
- `cutoff`
  - metal environment selection cutoff
- `additional_residue_ids`
  - residues to force into the model
- `charge_mode`
  - `resp` or `fixed`
- `qm_backend`
  - `pyscf` or `psi4`
- `basis`
- `scale_factor`
- `frcmod_files`
- `gaff`
- `force_field`
- `water_model`

## Proposed Result Object

`Xponge.MCPB()` should return a structured result, but the primary success
criterion is that the parent `Molecule` has been updated into a SPONGE-export-
ready state.

Suggested shape:

```python
result = {
    "molecule": updated_molecule,  # directly save_sponge_input-ready
    "small_model": small_model,
    "large_model": large_model,
    "frcmod_path": ".../metal_center.frcmod",
    "connect_records": [...],
    "updated_charge_atoms": [...],
    "registered_metal_templates": ["ZN"],
    "sponge_ready": True,
    "metadata": {
        "method": "seminario",
        "qm_backend": "pyscf",
        "charge_mode": "resp",
    },
}
```

## Mapping To Amber MCPB.py Concepts

Amber `MCPB.py` is step-based. The `Xponge` version should expose one workflow,
but internally it can still map to similar phases.

### Phase A. Environment selection

Amber-like concepts:

- `ion_ids`
- `additional_resids`
- `cut_off`
- `add_bonded_pairs`

Xponge draft:

- select target metal atoms
- select coordinating atoms and nearby residues
- optionally let the user override missing or ambiguous bonded pairs

### Phase B. Model construction

Amber-like concepts:

- small model
- large model
- frozen atoms
- charge and spin assignment

Xponge draft:

- build a compact QM-ready `small_model`
- build a larger context `large_model`
- adapt both to `QMMolecule`

### Phase C. Parameter generation

Amber-like concepts:

- empirical
- Seminario
- modified Seminario
- Z-matrix
- blank `frcmod`

Xponge draft first-pass methods:

- `blank`
- `empirical`
- `seminario`

Deferred:

- `modified_seminario`
- `zmatrix`

### Phase D. Charge fitting

Amber-like concepts:

- local RESP refit with different restraint levels
- fixed-charge fallback

Xponge draft:

- use the shared `qm/` backend for SCF and ESP
- use the existing C++ RESP core for numerical fitting
- update only the metal center and selected nearby ligand atoms/residues

### Phase E. Export integration

Amber-like concepts:

- mol2
- frcmod
- leaprc / final setup files

Xponge draft:

- write `frcmod`
- patch parent-molecule charges in memory
- add explicit metal bond/connect data to the parent `Molecule`
- ensure topology construction can proceed without requiring a prebuilt whole-
  complex residue template
- preserve or add metal `CONECT` on PDB export
- optionally write metal-centered local mol2 or debug artifacts

## Supported Modes In The Draft

### Mode 1. Bonded metal-center model

This is the main target for a MCPB-like workflow.

Expected outputs:

- a parent `Molecule` that can directly `Save_SPONGE_Input(...)`
- bonded `frcmod`
- local charge updates
- explicit connect records

### Mode 2. Nonbonded metal model

This should also be designed in from the start, even if implemented after the
bonded mode.

Expected outputs:

- nonbonded ion treatment
- optional local charge updates
- no bonded force-constant derivation

## Element Acceptance And Validation Tiers

### API acceptance

`Xponge.MCPB(...)` should accept any valid periodic-table element symbol in
`ion_info.element`.

This means:

- do not hard-code the public API to only a small list of metals
- allow all elements to enter request validation
- separate "accepted as input" from "validated for production use"

### Tier 1: First bonded-validation target

These are the best first candidates for validated bonded MCPB-like support:

- `Mg`
- `Ca`
- `Zn`
- `Fe`
- `Mn`
- `Co`
- `Ni`
- `Cu`

Reasons:

- common in biomolecular systems
- likely highest practical value
- more realistic first test set

### Tier 2: Broader likely bonded or nonbonded validation

- `Mo`
- `W`
- `V`
- `Cr`
- `Ru`
- `Rh`
- `Pd`
- `Ag`
- `Cd`
- `Pt`
- `Au`
- `Hg`

### Tier 3: Accept in API, but not guaranteed initially

- lanthanides
- actinides
- rare heavy-metal edge cases

For these, the API should still accept the element name, but MCPB should
require explicit validation of:

- available ion typing / template data
- bonded/nonbonded parameter availability
- QM basis or ECP availability
- charge / spin consistency

## File-format Risks To Handle Explicitly

### PDB element ambiguity

If the PDB element field is missing or malformed, current fallback rules may
mis-read metal names from atom names.

Implication:

- `MCPB` should prefer explicit user-provided `ion_info.element`
- PDB import should not be trusted as the sole source of metal identity
- API acceptance should not depend on PDB parser confidence

### MOL2 element/type ambiguity

Current `mol2` reading often derives element identity from atom type mass or
fallback name parsing.

Implication:

- metal atom type and element must be validated before MCPB continues
- if element resolution is ambiguous, stop early with a targeted error

### Assignment bond-order risk

Current whole-system `Assign` bond-order logic is not robust for arbitrary
metal centers.

Implication:

- do not start MCPB from a full-system `Assign`
- only use `Assign` for local, controlled fragments

## Proposed Internal Module Layout

Suggested new modules:

- `src/XpongeCPP/mcpb/__init__.py`
- `src/XpongeCPP/mcpb/api.py`
- `src/XpongeCPP/mcpb/models.py`
- `src/XpongeCPP/mcpb/selection.py`
- `src/XpongeCPP/mcpb/model_builder.py`
- `src/XpongeCPP/mcpb/charge_refit.py`
- `src/XpongeCPP/mcpb/frcmod.py`
- `src/XpongeCPP/mcpb/export.py`

Suggested responsibilities:

- `models.py`
  - MCPB request/result dataclasses
- `selection.py`
  - identify metal environment and enforce explicit checks
- `model_builder.py`
  - build small and large models
- `charge_refit.py`
  - local RESP/fixed-charge workflows via `XpongeCPP.qm`
- `frcmod.py`
  - generate bonded or placeholder frcmod content
- `export.py`
  - parent-molecule patching for SPONGE export readiness
  - connect records, charge patching, PDB/mol2 integration

## Proposed Export Behavior

### Primary export contract

After `Xponge.MCPB(...)` completes successfully in bonded mode:

- the returned `result["molecule"]` should be immediately usable with
  `Save_SPONGE_Input(...)`
- the same should also be true for the in-memory `molecule` object if MCPB is
  implemented as an in-place workflow
- users should not be required to manually assemble `frcmod + connect +
  charge-overrides` before SPONGE export

### Charge override behavior

After local RESP fitting:

- update the metal center atom charges
- update nearby ligand/residue charges
- preserve the rest of the system

### PDB export behavior

PDB writer should preserve or inject:

- metal-ligand `CONECT`
- correct metal element field

### Force-field export behavior

The workflow should generate a self-contained metal-center parameter artifact,
but this is secondary to making the live `Molecule` exportable:

- `frcmod`
- optional local mol2 or residue-template fragment
- metadata describing method and QM backend

### SPONGE export behavior

The MCPB workflow should patch the parent `Molecule` so that SPONGE topology
construction succeeds using:

- explicit metal-ligand bonds in `explicit_bonds` / `residue_links` as needed
- registered metal-ion template support for single-ion residues when applicable
- MCPB-generated bonded parameters available through the active force-field
  registries
- in-memory charge overrides already applied to atoms

## First Implementation Plan

### Phase 1. API and environment selection

- add `Xponge.MCPB(...)` placeholder API
- define request/result dataclasses
- support explicit `ion_ids`, `ion_info`, `bonded_pairs`
- perform strict validation on metal identity

### Phase 2. Model extraction

- extract metal-centered residue neighborhood from `Molecule`
- build small and large model objects
- add tests on a simple Zn-binding example

### Phase 3. Charge refit

- local RESP route using shared `qm/`
- charge patching back into the parent `Molecule`

### Phase 4. Parent-molecule integration

- register or validate metal-ion template support needed by export
- add explicit metal-ligand connect/bond information to the parent `Molecule`
- ensure the parent `Molecule` becomes `Save_SPONGE_Input(...)`-ready

### Phase 5. frcmod generation

- add `blank`
- add `empirical`
- add `seminario`

### Phase 6. Export integration and artifacts

- PDB `CONECT`
- charge override persistence
- artifact writing
- direct `Save_SPONGE_Input(...)` validation on MCPB-treated systems

## Suggested First-Release Scope

To keep the project moving, the first useful `Xponge.MCPB()` should probably be:

- bonded model only
- all elements accepted at the API layer
- Tier-1 metals validated first for bonded production workflows
- explicit `ion_ids`, `ion_info`, and `bonded_pairs` required
- RESP-based local charge refit
- `blank` + `seminario` force-constant modes
- a returned `Molecule` that can directly `Save_SPONGE_Input(...)`

This is narrower than full Amber MCPB.py parity, but it would already cover the
most valuable real cases and fit the current `XpongeCPP` architecture cleanly.

## Open Questions

- Should `ion_ids` refer to global atom ids only, or also support residue-local
  selectors?
- Should `bonded_pairs` be mandatory in v1 for bonded mode?
- Should local RESP modify whole ligating residues or only the selected atoms?
- Should MCPB mutate the input `molecule` in place, or return a patched copy by
  default?
- Should the first export target stop at direct SPONGE export readiness, or
  also provide a one-shot PDB-plus-artifact packaging helper?
- For non-Tier-1 elements, should unsupported basis/ECP combinations fail early
  in `qm` validation or inside `MCPB`?
