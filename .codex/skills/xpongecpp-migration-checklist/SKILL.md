---
name: xpongecpp-migration-checklist
description: Use when planning or executing migration of an Xponge workflow into XpongeCPP, especially Amber protein-ligand-solvent workflows such as 1KV2+B96. Provides the parity checklist, baseline procedure, implementation order, and acceptance criteria.
---

# XpongeCPP Migration Checklist

## Scope Rule

Migrate workflows, not Xponge's Python object internals. Preserve common Python API names only where needed for user scripts.

The reference implementation is the original repository at:

`/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge`

Treat it as read-only. All implementation, tests, and docs belong in:

`/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/XpongeCPP`

## Baseline Procedure

For every migrated workflow:

1. Create a minimal Xponge script that runs the original workflow.
2. Use the same input files, random seeds, solvent distance, ion counts, and parameter files.
3. Save SPONGE input and any PDB/mol2 output to `/tmp/xponge_*`.
4. Create the equivalent XpongeCPP script and save to `/tmp/xpongecpp_*`.
5. Compare counts and numeric records.
6. Add tests that either:
   - generate the local Xponge baseline when available, or
   - compare against locked expected headers/counts and canonical numeric records.

Use the correct original Xponge parameter registration API. For Amber `frcmod`, this is usually:

```python
from Xponge.forcefield import amber
amber.load_parameters_from_frcmod(path, prefix=False)
```

Do not confuse this with top-level `Xponge.load_frcmod`, which parses data but may not register parameters in the force-field registry.

## Migration Completion Checklist

- Python import path exists and registers templates and parameters.
- Parser preserves all source information needed by topology.
- C++ data model stores explicit connectivity and residue links.
- Topology generation follows Xponge's original algorithm.
- Exported SPONGE headers match Xponge.
- Numeric records match Xponge or have a documented equivalence transform.
- Safety tests cover malformed inputs and out-of-range IDs.
- Performance-critical loops are in C++.
- Full test suite passes.

## Amber Protein-Ligand Workflow Checklist

For workflows such as `1KV2_H.pdb + B96.mol2 + B96.frcmod + ff14SB + GAFF + TIP3P + ions`:

- Import `XpongeCPP.forcefield.amber.ff14sb`.
- Import `XpongeCPP.forcefield.amber.gaff` or `gaff2`.
- Import the solvent model module.
- Register ligand `frcmod`.
- Load protein PDB.
- Load typed ligand mol2.
- Merge or assemble protein + ligand without losing mol2 connectivity.
- Add solvent box.
- Replace selected water residues with ions.
- Save PDB and SPONGE input.
- Compare with Xponge baseline.

## Required Comparisons

Always compare:

- atom count and residue count;
- residue names and atom names;
- atom type names;
- coordinates and box;
- mass and charge;
- LJ type mapping and coefficients;
- bond count and numeric records;
- angle count and numeric records;
- proper and improper dihedral records;
- exclude list;
- nb14 list.

Canonicalize only where Xponge's force class defines same-force equivalence.

## Do Not Migrate Yet Unless Requested

- Full Xponge dynamic Python object graph.
- Non-Amber force fields.
- GUI-specific helpers.
- MindSponge runtime integration.
- FEP/special force-field modules.
- Byte-for-byte output formatting beyond what the workflow requires.

## Stop Conditions

Stop and investigate before coding if:

- Xponge and XpongeCPP disagree but no Xponge source path has been identified.
- A proposed fix adds chemistry heuristics not present in Xponge.
- A fix only changes output sorting without proving the order is semantically irrelevant.
- A missing parameter would be hidden by a default guessed parameter.
- A molecule rebuild drops explicit bonds, residue links, or box state.
