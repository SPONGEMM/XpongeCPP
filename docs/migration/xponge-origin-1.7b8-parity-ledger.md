# Xponge-origin 1.7b8 Parity Ledger

## Baseline

| Repository | Revision | Role |
| --- | --- | --- |
| Xponge-origin | `1.7b8` / `695ea0b4b8555d545cdb935545b076b018bdddd7` | Compatibility target |
| XpongeCPP | `0.2.0` development tree | Unreleased implementation |

The earlier 1.7b7 ledger remains the historical record for bundle, RESP,
Hessian, empirical fitting, force-field data, and MCPB removal. This ledger
records the incremental 1.7b8 behavior change.

## 1.7b8 delta

| Origin behavior | XpongeCPP implementation | Status |
| --- | --- | --- |
| Hash-closed `MetalParameterPatch` wire contract | `metal_assignment/patch.py` is source-compatible with origin | implemented |
| Derived-model-only `MetalLocalModelPackage` | input/provider modules synchronized with origin | implemented |
| Existing ordinary molecule owns the base force field | bonded fit supports `base=None` and `preassigned_molecule` provenance | implemented |
| `ParameterizationResult` never represents an applied complete system | result statuses, completeness, and audit validation synchronized | implemented |
| `apply(molecule, patch, ...)` changes only local parameters | C++ atom state, molecule-local LJ/bond/angle overrides, and unique-type dihedral registrations | implemented |
| Transactional registry and molecule rollback | Python and native force-field registries are snapshotted and restored on failure | implemented |
| Embedded metal residue preparation | dynamic residue templates and `prepare_residue_templates` supported | implemented |
| Assigned-system and apply-input APIs removed | modules, exports, and legacy public signature removed | implemented |
| Nonbonded full-system assignment rejected by molecule API | origin rejection test synchronized | implemented |
| PDB chain changes split equal residue numbers without `TER` | already native in the C++ reader; origin regression test added | implemented |
| Missing template head/tail atom does not abort chain inference | already guarded by C++ atom lookup | implemented |

## Verification gates

- Origin 1.7b8 contracts, artifacts, and local patch application tests.
- RESP, Hessian, base assignment, standard assignment, and provider tests.
- Local-patch raw/bundle/legacy real-SPONGE numerical gates for constrained
  RESP, manual bonded fitting, and Seminario fitting.
- XpongeCPP API regression for equal residue numbers across PDB chains.
- Release metadata and namespace checks for unreleased version `0.2.0`.

The old assigned-system numerical fixture and its public helpers are not a
compatibility surface. They were replaced by a local-patch fixture rather than
retained as deprecated APIs.
