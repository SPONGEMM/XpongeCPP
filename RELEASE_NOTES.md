# XpongeCPP 0.2.16

Compatibility target: Xponge-origin 1.7.7.

Add public `get_explicit_topology`, `set_explicit_topology`, and `set_atom_lj`
interfaces for validated external parameters, including precomputed ByteFF
topologies. Preserve ordered periodic rows, explicit 1-4 scales and exclusions;
use content-addressed LJ types to isolate molecule parameters. The native FEP
interpolation algorithm is unchanged in this release.

Installed-wheel release tests now cover the API, materialized topology copies,
disjoint assembly, full/partial MCS projection, and FEP energy interpolation.
See `docs/explicit-topology.md` for units, subset scope and limitations.
No ByteFF model inference, real FEP or full-engine force validation is claimed.

# XpongeCPP 0.2.15

Compatibility target: Xponge-origin 1.7.6.

Includes the public metal-patch Atom endpoint fix, deterministic AMBER improper
ordering and retained dummy bonded scaffolds described below. Complete the
synthetic CMAP and water-link test parameters and use NALA + internal ALA for
the cross-component peptide fixture. Missing production parameters still fail
strictly; no default force-field terms are introduced.

The CDK8 parameter-provenance and cross-backend ordered-improper limitations
below remain; full-engine/MD equivalence is not claimed.

# XpongeCPP 0.2.14 (unpublished candidate)

Compatibility target: Xponge-origin 1.7.5.

Metal patch link validation now accepts public ResidueLink Atom endpoints as
well as integer atom indices. Missing required links are still rejected; a
regression test covers both cases.

## deterministic AMBER impropers (2026-10-23)

Native periodic improper matching retains original parameter slots and all
selected periodicities, with latest per-period override and stable
type/residue/local-name ambiguity resolution. Parameter lookup caches retain
slot metadata rather than occurrence-dependent permutations. Materialized FEP
endpoints continue to bypass rematching. See AMBER_IMPROPER_POLICY.md; this
does not emulate LEaP's ambiguous unit-cache ordering or change harmonic
improper definitions. No frozen input or trajectory is automatically migrated.

## retained dummy bonded scaffolds (2026-10-23)

Dual construction retains opposite-branch source bond, angle and ordered
periodic torsion entities touching non-MCS atoms at both endpoints. Ghosts still
have zero LJ/charge, and endpoint NB14 uses these nonbonded states. Common-only
terms keep their physical endpoint parameters. No union-graph reparameterization
or extra cross-branch angles/torsions is performed. Ordinary modeling and the
explicit low-level one-sided soft-bond interface are unchanged.

This corrects the physical-only projection described below: auxiliary dummy
terms are part of the chosen dual-topology protocol, not unwanted physical
endpoint contamination. It is not a globally validated replacement. CDK2
historical ligand-subset tables match at seven step2 states; CDK8 still has
historical parameter-provenance and cross-backend ordered-improper discrepancies.
Do not reuse old trajectory acceptance or claim full engine/MD equivalence.

# XpongeCPP 0.2.13

Compatibility target: Xponge-origin 1.7.4.

Disjoint molecule assembly now preserves materialized physical topologies from
either operand and appends all offset force/exclusion indices, including newly
added solvent. Ordinary uncached assembly remains lazy. A requested covalent
link between cached topologies without an explicit combined physical topology
is rejected before mutation rather than silently retaining stale force tables.

Soft-core legacy export omits the empty ordinary NB14 table, whose native
loader requires ordinary LJ tables that soft-core mode intentionally omits.
Explicit NB14-extra coefficients are retained. Nonempty ordinary NB14 terms
are exported as explicit coefficients using exactly their ordinary A-state LJ
parameters and scales, preserving legacy manual softcore-helper behavior without
mutating the input molecule. Already interpolated FEP extras are not reweighted;
ordinary non-FEP LJ/NB14 export is unchanged.

FEP dual-state construction now projects materialized physical endpoint force
entities and exclusions into union atom IDs, preserving ordered improper terms
instead of selecting them again from a union bond graph. The unchanged
environment is projected separately for the B state. Covalently linked target
or replacement residues without an explicit physical B environment topology
are rejected rather than assigned an invented endpoint. This safety boundary
only affects FEP dual construction; ordinary modeling is unchanged.

Materialized FEP export now groups ordinary hard-bond components within each
residue, returning the full atom-identity permutation. Ghost/soft bonds are not
used to bypass the SPONGE continuity contract. Pure atom/residue permutations
preserve and remap cached bonds, angles, ordered periodic torsions, exclusions
and NB14 terms, rather than dropping the materialized topology. Non-FEP default
atom ordering and the continuity rejection guard are unchanged.

FEP bonds present in only one physical endpoint now occur only in the soft-bond
table, matching Xponge's bond merge partition. They are no longer duplicated as
lambda-weighted ordinary harmonic bonds. Soft parameters come from the actual
endpoint topology, including overrides, rather than a parameter/distance lookup.
Pairs present in both endpoints with different equilibrium lengths retain their
ordinary lambda-weighted terms. Ordinary non-FEP topology/export is unchanged.

Amber parmdat registration now stores the atom-section masses in the native
parameter registry, matching Xponge rather than silently relying on a partial
element fallback. MOL2 loading and template registration use the supplied mass;
later parmdat registrations and explicit frcmod MASS overrides retain their
normal precedence. No element-specific output correction is applied.

PDB and mmCIF loading now gives file-declared inter-residue connections
precedence over an additional implicit template link between two connection
sites that already have explicit partners. This prevents branched multi-residue
structures from acquiring unintended head/tail cycles. The rule does not use
residue names, coordinate distances, or ATOM/HETATM classification. Implicit
polymer and multi-residue lipid links remain enabled when the two connection
sites are not both explicitly assigned; explicitly declared neighbor links and
ring closures are retained. PDB `ignore_conect=True` still ignores CONECT.

Xponge-origin 1.7.4 uses the same explicit-link precedence. mmCIF edit operations
still run after inference, and external residue-link overrides keep their
final-overlay semantics. Missing connectivity still requires explicit records
or suitable templates.

mmCIF loading now accepts `ignore_conect=True` and filters noncovalent
`struct_conn`, Mokda semantic, and bond-edit records. The ignore option bypasses
file bond declarations and edits while retaining template bonds and separately
supplied residue links. PDB and mmCIF therefore support the same Mokda-controlled
connection policy without requiring a complete-graph interface.

The read-only `amber_parameter_exists` query checks registered Amber parameter
coverage without mutating the registry. Missing bond and angle parameters now
raise an error instead of generating arbitrary defaults. Bare TIP3P registration
supplies its actual model bond parameters and zero-strength angle entries.

# XpongeCPP 0.2.12

Add official Amber OL24 DNA parameters and reproducible OFF-derived residue
templates, including HOP3-bearing 5-prime monophosphate templates. OL3 and OL24
register their own RNA/DNA terminal templates without an additional user import.
Preserve default hydroxyl termini and document parameter provenance and citations.

Fix nucleic polymer anchors and water-specific atomic ion aliases. This release
also includes full mmCIF chain identity and contiguous bundle component ordering
with source-atom and protocol remapping.

# XpongeCPP 0.2.11

Preserve full mmCIF chain identifiers in residue grouping, chain-end inference,
explicit terminal selectors and external residue links. Distinct chains such as
ML1 and ML2 no longer collapse to M or acquire unintended cross-chain bonds.
Python residue identities and template replacement retain the full chain name.
PDB reading and writing behavior is unchanged.

# XpongeCPP 0.2.10

Fix bundled SPONGE input export to group linked residue components in contiguous
atom order, matching the existing raw-text exporter. Both Python and native
bundle entry points validate molecular component continuity before writing topology.

Source atom mappings, collective variables (including typed virtual atoms),
distance constraints, positional reference coordinates and SITS selections follow
the reordered atoms. Selection weights and named virtual-atom references retain
their original meaning. Regression tests cover raw/native/Python bundle export,
covalent and coordination links, source mappings and protocol remapping.

# XpongeCPP 0.2.8

Compatibility target: Xponge-origin 1.7.1.

This patch release restores the Python instance-attribute behavior required
by the Xponge RESP diagnostics contract. The pybind-backed `Assign` type now
supports compatibility attributes, allowing `calculate_charge(...,
return_diagnostics=True)` and `return_metadata=True` to retain
`charge_fit_diagnostics` and `charge_fit_metadata` after fitting.

Regression coverage verifies diagnostics, metadata, and general legacy
extension attributes without requiring a QM backend. Release validation also
includes a real PySCF RESP calculation and the complete compatibility suite
against the official Xponge 1.7.1 package.

# XpongeCPP 0.2.7

Compatibility target: Xponge-origin 1.7b11.

This patch release adds typed virtual-atom support to the bundled SPONGE
protocol path and restores legacy SPONGE DAT handling in the MDAnalysis
compatibility layer. It also includes the release-CI test-path correction so
the Xponge-origin compatibility cohort is collected and executed during
package validation.

# XpongeCPP 0.2.6

Compatibility target: Xponge-origin 1.7b10.

This parity release changes QM backend selection from an operating-system rule
to runtime capability detection. XpongeCPP now prefers an installed PySCF
backend, falls back to Psi4 when PySCF is absent, and reports a clear error
when neither backend is available. Windows packaging continues to skip PySCF;
users may supply a compatible PySCF installation or install Psi4 separately.

Regression coverage verifies PySCF preference, Psi4 fallback, and the
no-backend error path. Windows wheel CI remains a packaging smoke test; it does
not claim native Windows PySCF support.

# XpongeCPP 0.2.5

Compatibility target: Xponge-origin 1.7b9.

This patch preserves source-atom identity when raw SPONGE export reorders
noncontiguous linked residue components. The native layer now reports the
saved atom permutation before serialization, and the compatibility layer
applies the same permutation to caller-provided source IDs. This keeps Mokda
trajectory topology and atom-order mappings aligned for proteins, ligands,
metals, solvent, and repeated solvent molecules without matching atoms by
name or coordinates.

Focused release validation covers covalent and coordination-driven residue
reordering, source-ID mapping, and the complete native bundle regression
cohort. The complete repository suite passes with 716 tests, 10 optional
skips, and 1 expected xfail. The change adds one linear residue/atom
permutation pass before raw save and does not move force-field assignment back
to Python.

# XpongeCPP 0.2.4

Compatibility target: Xponge-origin 1.7b9.

This compatibility patch preserves residue-link bonds in public MOL2 exports,
keeps mixed standard/GAFF systems intact when invoking GAFF or GAFF2
parameterization, and aligns the legacy MDAnalysis raw-SPONGE reader with the
format hints exposed by Xponge. These fixes cover the Mokda SPONGE-input path
without moving force-field assignment back to Python.

Release validation includes the focused MOL2, mixed-force-field, and
MDAnalysis compatibility regressions. The complete repository suite passes
with 715 tests, 10 optional skips, and 1 expected xfail. A CPython 3.11 wheel
was built and imported from an isolated environment under the Mokda resource
gate with `OOMKillDelta=0`.

# XpongeCPP 0.2.3

Compatibility target: Xponge-origin 1.7b9.

This release completes the Mokda-used same-script compatibility surface and
adds scientific-output comparison for real systems. It preserves native
template bonds when an mmCIF residue also supplies explicit covalent bonds,
including model pseudo-bonds such as the TIP3P H-H distance constraint, while
still honoring explicit Mokda edit deletions. It also synchronizes temporary
legacy residue-link overrides with the native registry and canonicalizes CMAP
grid identities independently of force-field declaration order. The legacy
MDAnalysis raw-SPONGE reader now advertises the `SPONGE_MASS` format and
recognizes `_mass.txt`, matching the Xponge interface used by Mokda H5MD
analysis.

Linux x86_64 release-candidate evidence:

- XpongeCPP repository suite under the coherent managed 0.2.3 runtime after the
  raw-reader compatibility fix: 717 passed, 7 optional skips, 1 expected
  xfail;
- Xponge-origin 1.7b9 tests through the installed compatibility namespace:
  36/36 passed;
- Mokda auto-discovered Xponge consumer matrix: 327 passed, 14 conditional
  skips, and 12 subtests under both providers;
- real 3GOU/4EWL Life Modeling under the isolated XpongeCPP 0.2.3 managed
  runtime completed in 205.656 s with a 1.09 GiB peak and
  `OOMKillDelta=0`;
- complete Mokda Lipid21/1BL8 tier: 8/8 files passed, with a 1.995 GiB peak
  and `OOMKillDelta=0`;
- provider-neutral scientific manifests match for repaired 1EMA, 1KV2+B96,
  8RYK, 1BNA, 1GYA, and 1BL8; the existing 3GOU/4EWL four-quadrant metal
  acceptance also passes under both managed providers;
- provider-neutral 1KV2 benchmark remains 42,600 atoms in both
  implementations; XpongeCPP median 0.0588 s versus Xponge 12.47 s, and the
  measured XpongeCPP pre/post change is +2.57%, within the 5% release
  tolerance.

# XpongeCPP 0.2.2

Compatibility target: Xponge-origin 1.7b9.

This compatibility release closes the Mokda-used Xponge surface under the
same-script `import Xponge` contract. It synchronizes legacy atom-type and
residue-link state with the native registries, aligns GAFF/GAFF2, PDB/link,
Lipid17/Lipid21 and non-Amber behavior with the current reference, and ships
the non-Amber reference force-field data inside wheels. The release also adds
a fail-closed machine-readable capability manifest and validates both source
and installed-wheel execution.

Linux x86_64 release-candidate evidence:

- installed CPython 3.12 wheel: 701 passed, 7 optional skips, 1 expected
  unsupported legacy-frcmod workflow;
- Xponge-origin 1.7b9 suite through the compatibility namespace: 36/36 passed;
- Mokda auto-discovered Xponge consumer matrix: 324 passed, 13 conditional
  skips and 12 subtests;
- provider-neutral 1KV2 benchmark: 42,600 atoms from both implementations;
  XpongeCPP median 0.0588 s versus Xponge 12.47 s;
- pre/post XpongeCPP micro-benchmark change: +2.57% total, within the 5%
  release tolerance; no assignment hot path was moved back to Python.

# XpongeCPP 0.2.1

Compatibility target: Xponge-origin 1.7b9.

This patch release corrects GLYCAM terminal-zero residue metadata across the
native template registry, bundled Python data, and the pinned reference
force-field copy. Terminal-zero names such as `0MA`, `0aA`, `0AD`, and `0aD`
now carry no synthetic `O0`/`C0` head attachment. The release adds parity and
standalone PDB-export regression coverage for all four pyranose/furanose and
D/L representative families.

# XpongeCPP 0.2.0

Compatibility target: Xponge-origin 1.7b8.

This is a breaking feature-alignment release from XpongeCPP 0.1.6.
Supported Python versions are CPython 3.10 through 3.12.

## Highlights

- Ports the SPONGE bundle contract, native saver, legacy/bundle converters,
  restart and output handling, H5MD readers, MDAnalysis integration, and CLI
  entry points introduced through Xponge-origin 1.7b8.
- Exposes generic RESP linear constraints, equivalence groups, constraint
  preflight, RHF/ROHF/UHF reference selection, SCF strategy controls, progress
  events, and fit metadata through the public charge API.
- Replaces the former MCPB-specific implementation with the origin-style
  `XpongeCPP.metal_assignment` package. The package includes immutable
  contracts and artifacts, base and nonbonded assignment, constrained RESP,
  Hessian providers, empirical/manual/Seminario bonded fitting, hash-closed
  local parameter patches, transactional application to preassigned
  molecules, provenance, and SPONGE export mappings.
- Aligns application semantics with Xponge-origin 1.7b8: removes the temporary
  assigned-system/apply-input materialization layer and applies
  `MetalParameterPatch` objects without replacing the ordinary force field.
- Removes the `XpongeCPP.MCPB` API and the `XpongeCPP.mcpb` package. Callers
  must migrate to `XpongeCPP.metal_assignment`.

## Verification snapshot

On Linux x86_64 with CPython 3.11:

- bundle contract/conversion/analysis cohort: 166 passed with the discovered
  SPONGE executable;
- after adding environment, `PATH`, and sibling-build discovery, the focused
  bundle converter/SPONGE zero-frame gate passed all 7 tests with the explicit
  environment override unset;
- rebuilt RESP, bundle, and metal cohort: 238 passed, 12 optional QM/backend
  tests skipped;
- metal-assignment origin-parity cohort: 82 passed, 23 optional
  Mokda-chemcore/backend tests skipped;
- final combined pre-1.7b8 release/RESP/bundle/metal/numerical parity cohort:
  217 passed, 23 optional external-provider tests skipped;
- Xponge-origin 1.7b8 patch/contract/provider/release cohort:
  85 passed, 16 optional external-provider tests skipped;
- RESP constraints, charge models, metal RESP, base-charge, and Hessian cohort:
  99 passed, 11 optional backend tests skipped;
- 1.7b8 local-patch constrained RESP, explicit-reference manual bonded, and
  deterministic-Hessian Seminario real-SPONGE gates pass raw-to-bundle-to-raw
  numerical comparison;
- wheel and source distribution built successfully and passed `twine check`;
- a clean wheel smoke environment imported `XpongeCPP`,
  `XpongeCPP.io_bundle`, and `XpongeCPP.metal_assignment`, reported version
  `0.2.0`, and exposed no `MCPB` symbol.

The repository-wide suite currently has 504 passing tests and 31 known
unrelated failures, compared with 412 passing and 32 failing in the pre-sync
baseline. The failures are tracked separately and do not occur in the bundle,
RESP, or metal-assignment parity cohorts. An actual cross-platform CI run and
a bundle-capable SPONGE runtime remain release gates rather than claims made by
this local verification snapshot. The latter has an explicit
`SPONGE_BUNDLE_EXECUTABLE`-controlled acceptance test and is not silently
treated as passing on legacy-only runtimes.
