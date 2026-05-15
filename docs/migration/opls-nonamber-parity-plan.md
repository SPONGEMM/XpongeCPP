# XpongeCPP OPLS/Non-Amber Parity Migration Plan

## Summary

This plan aligns `XpongeCPP` with the current `Xponge` behavior for the OPLS/non-amber
workflow and closes three previously incomplete historical compatibility gaps:

1. OPLS/non-amber EMC polymer workflow parity using the `poly` regression case
2. `19b2c4b`: ITP loading when the same residue name is created multiple times
3. `e28bd6b`: `save_mol2` bond mapping by atom identity rather than atom name
4. `032d574`: hybrid-36 structure saving for molecule and assignment PDB writers

The acceptance standard is strict file-level parity with `Xponge` whenever `Xponge`
already defines a stable output.

## Key Changes

### 1. OPLS initialization and parameter registration

- Replace the current `XpongeCPP.forcefield.opls` placeholder namespace with a true
  OPLS compatibility initializer.
- Mirror the behavior of `Xponge.forcefield.opls.__init__`:
  - OPLS LJ combining rules
  - OPLS bond/angle/improper/RB type-name getter semantics
  - 1-4 exclusion policy
  - `load_parameter_from_ffitp(filename, folder, reset=True)` compatibility entrypoint
- Expose `load_parameter_from_ffitp` from the top-level `XpongeCPP` namespace.

### 2. `poly` non-amber workflow parity

- Use the vendored test data under `tests/data/poly` as the stable regression fixture.
- Make the following workflow match `Xponge` output file-for-file:
  - `load_parameter_from_ffitp(ps.itp.name, ps.itp.parent)`
  - `load_molpsf(ps.psf)`
  - normalize the trailing 9-field GRO box line to 3 fields
  - `load_gro(ps.gro, system)`
  - `save_sponge_input(system, prefix)`
- Match all generated SPONGE files, not only topology structure.

### 3. `19b2c4b` repeated residue compatibility

- Port `Xponge`’s residue-type reuse logic for repeated residue names in ITP loading.
- When a residue name reappears with an atom type or charge mismatch:
  - first try to reuse an already-created compatible residue type
  - only create a new derived residue type if no compatible one exists
- Preserve the original `Xponge` distinction between head/tail cases and normal cases.

### 4. `e28bd6b` `save_mol2` atom-id mapping

- Ensure `save_mol2` maps residue connectivity through the actual residue atom instance,
  not through atom names.
- The implementation must remain correct when atom names are duplicated or ambiguous
  inside one residue.

### 5. `032d574` hybrid-36 saving

- Add hybrid-36 integer formatting for:
  - molecule PDB atom serials
  - molecule PDB residue sequence numbers
  - LINK/SSBOND/CONECT fields
  - assignment PDB atom serials and CONECT fields
- Match `Xponge` overflow behavior when the hybrid-36 range is exceeded.

## Implementation Order

1. Implement the OPLS initialization layer and `load_parameter_from_ffitp`.
2. Fix `poly` workflow parity and remove its current `xfail`.
3. Port the `19b2c4b` repeated-residue ITP logic.
4. Fix `save_mol2` atom identity mapping for `e28bd6b`.
5. Add hybrid-36 saving parity for `032d574`.
6. Re-run the relevant non-amber, MOL2, and PDB regression suites.

## Test Plan

### `poly` workflow

- `tests/test_poly_non_amber_workflow.py`
- Remove the current `xfail` once the generated outputs match `Xponge`.

### `19b2c4b`

- Add a dedicated ITP regression with the same residue name reused under conflicting
  atom typing or charge.
- Verify that the resulting residue-type selection matches `Xponge`.

### `e28bd6b`

- Add a dedicated `save_mol2` regression using duplicated or ambiguous atom names in a
  residue.
- Verify that connectivity in the exported MOL2 matches the original topology and the
  `Xponge` reference output.

### `032d574`

- Add PDB saving tests that explicitly validate:
  - atom serial hybrid-36 transitions
  - residue sequence hybrid-36 transitions
  - LINK/CONECT hybrid-36 fields
  - overflow errors
  - assignment PDB writer parity

### Regression protection

- Re-run the relevant existing tests:
  - `tests/test_non_amber_parsers.py`
  - `tests/test_psf_io.py`
  - `tests/test_b96_mol2_gaff.py`
  - `tests/test_assign_charge_models.py`
  - `tests/test_pdb_chain_terminal_semantics.py`

## Assumptions and Defaults

- `Xponge` remains the single behavior reference for this migration.
- Parity means strict textual equality for stable exported SPONGE files unless a test
  explicitly encodes a looser canonical comparison.
- If a fix belongs in the C++ core instead of the Python compatibility layer, fix it in
  the core rather than adding Python-side output rewriting.
