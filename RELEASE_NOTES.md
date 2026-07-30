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
