# AMBER periodic improper policy: amber-deterministic-v1

The native AMBER topology builder uses the same convention as Xponge-origin:

1. Keep the central atom in parameter slot 3 (zero-based slot 2).
2. Match peripheral permutations against original parameter slots; prefer the
   fewest wildcards. At that specificity, the latest definition of each
   periodicity wins, and distinct periodicities remain separate terms.
3. Resolve ambiguous matches by the lexicographically least sequence of
   `(actual atom type, residue position, unique local atom name)` in parameter
   slot order, not global AtomId or candidate traversal order. Reject duplicate
   local names in a matching group rather than silently choosing an ID.
4. Cache original parameter slots and provenance order, not a previously chosen
   atom permutation. Repeated groups and exporter renumbering must not alter
   the function. Stable names/template identities are required; arbitrary
   renaming is outside this guarantee.

This is a deterministic AMBER-compatible convention, not exact reproduction of
LEaP's occurrence/cache-dependent sOrder behavior. The convention applies to
AMBER periodic improper construction, not CHARMM harmonic improper semantics.
The OPLS loader explicitly registers its shared periodic parameters with
`deterministic=False`, preserving its previous matching and ordering path.
Mixing AMBER and legacy conventions for the same matching group is rejected.
Reference discussion: https://github.com/openmm/openmm/issues/220

FEP builds and projects physical endpoints before mixing; it must not reselect
impropers on union IDs. Ordered terms, periodicities and phases are retained so
that U(lambda) = (1-lambda)*U_A + lambda*U_B. Only identical functions (including
whole reversal) are coalesced. Other equivalent permutations may remain as
separate weighted functions. Materialized topology overrides bypass matching.

Validate on nonplanar energy/gradient examples, several lambda values, repeated
groups, atom reorderings and cross-backend serialized identities. These are
local modeling checks, not full-engine or MD/free-energy acceptance. No frozen
input/trajectory is rewritten or automatically rerun by this source change.
