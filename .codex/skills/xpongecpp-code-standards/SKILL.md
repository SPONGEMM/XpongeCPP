---
name: xpongecpp-code-standards
description: Use when modifying XpongeCPP C++ core, Python bindings, topology, force-field parsing, solvation, assign, or SPONGE export code. Enforces Xponge parity-first implementation, flat C++ data structures, strict safety checks, and no ad-hoc topology rules.
---

# XpongeCPP Code Standards

## Non-Negotiable Rule

For topology, force matching, exclusions, nb14, solvation placement, ion replacement, and file export, preserve Xponge behavior unless the task explicitly says otherwise.

Do not invent extra chemical heuristics or special-case rules to make one test pass. If output differs from Xponge, first identify which original Xponge function produced the behavior, then port that algorithm or document a deliberate incompatibility.

## Required Workflow

1. Read the relevant Xponge implementation in `/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Xponge`.
2. Write or update a failing parity test before changing production code.
3. Implement the smallest C++ change that matches the Xponge behavior.
4. Compare parsed numeric output, not only text, unless byte-level compatibility is the explicit target.
5. Run the targeted test, then the full test suite.
6. Commit only coherent, verified changes.

## Architecture Boundaries

- `cpp/core/`: flat data model, ID validation, `Molecule`, `Residue`, `ResidueType`, `Assign`.
- `cpp/io/`: parsers and exporters only. Do not build force-field logic here.
- `cpp/forcefield/`: parameter registries, Amber parsing, template registration.
- `cpp/topology/`: bond graph, angle/dihedral/improper/exclude/nb14 generation.
- `cpp/solvation/`: solvent placement and ion/water replacement.
- `cpp/python/`: pybind11 bindings and stable API surface.
- `python/XpongeCPP/`: compatibility aliases and force-field import modules.

Keep hot loops in C++. Python wrappers should not iterate over atoms, bonds, waters, or topology terms.

## Data Model Rules

- Use `AtomId` and `ResidueId`; do not store raw cross-structure pointers.
- Keep `Molecule.atoms` and `Molecule.residues` contiguous.
- Keep topology independent from atom/residue storage.
- Preserve `Molecule.explicit_bonds` whenever rebuilding molecules.
- Validate bounds before public access; hot paths may use unchecked indexing only after prevalidation.
- Any batch mutation must leave `Molecule::validate()` true.

## Topology Parity Rules

- Follow Xponge's bonded-force construction pattern:
  - build connectivity/link distances;
  - generate candidate bonded terms from `topology_like` and `topology_matrix`;
  - group terms by `Same_Force`;
  - search exact types first;
  - search wildcard `X` patterns using Xponge's order and least-wildcard behavior.
- Proper dihedral multi-term reset semantics must match Amber/Xponge parsing.
- Improper dihedral must follow Xponge's `1-3-2-3` topology and same-force permutation around atom index 2.
- Exclude and nb14 must come from the same graph semantics as Xponge.
- Do not use distance guessing when explicit mol2 connectivity or residue template connectivity exists.
- Missing parameters should be treated as errors for parity work unless the test explicitly covers a Xponge fallback.

## Parser and Force-Field Rules

- `load_mol2` and template registration must share one parser path.
- Preserve mol2 residue names, atom names, atom types, charges, coordinates, and all declared bonds.
- `frcmod`/`parmdat` parsing must preserve Xponge/Amber ordering and override semantics.
- Public force-field imports in `python/XpongeCPP/forcefield/amber/` must register all required templates and parameters for that force-field path.

## Testing Requirements

- Add a minimal unit test for the exact behavior being changed.
- For parity tasks, generate or reuse a local Xponge baseline and compare:
  - headers/counts;
  - atom/residue names;
  - charges/masses/LJ;
  - bond/angle/dihedral/exclude/nb14 numeric meaning.
- Use canonical equivalence only when Xponge itself treats permutations as the same force.
- Always run:

```bash
rtk .venv/bin/python -m pytest -q
```

## Common Red Flags

- New `if residue.name == ...` branches outside clearly scoped compatibility code.
- New hard-coded bond, angle, or dihedral constants.
- Distance-based bond guessing used for typed mol2 or templated residues.
- Rebuilding `Molecule` without preserving `explicit_bonds` and `residue_links`.
- Python loops over large atom/water/topology arrays.
- Passing tests by sorting away semantically meaningful order without proving Xponge order is irrelevant.
