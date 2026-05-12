# Xponge vs XpongeCPP Architecture and Migration Status

This document compares the original Python Xponge package with the current
XpongeCPP repository. It is intended as a living migration map: it records what
is supported, what is partially supported, and what is still missing.

Status labels:

- **Supported**: implemented in XpongeCPP and covered by local tests or current
  regression fixtures.
- **Partial**: implemented for common or fixture-backed paths, but not yet
  equivalent to the full original Xponge behavior.
- **Not supported**: no usable XpongeCPP equivalent exists yet.
- **Not planned for core**: intentionally outside the current C++ core target.

## 1. Original Xponge Package Structure

The original Xponge package is centered on a small number of core modules, but
several support packages are also part of the modeling semantics.

| Xponge path | Role in original Xponge | Migration priority |
| --- | --- | --- |
| `__init__.py` | Public namespace, default constants, global aliases, default SPONGE writer registration. | High |
| `helper/` | Dynamic object system: `AtomType`, `ResidueType`, `Atom`, `Residue`, `ResidueLink`, `Molecule`, `GlobalSetting`, force type factories, `Xdict`, file/math helpers. | High |
| `assign/` | Small molecule assignment, bond order, ring markers, GAFF/SYBYL rules, TPACM4, RESP, Gasteiger, pH model, CIF/PDB/MOL2/XYZ/SMILES/PubChem entry points. | High |
| `forcefield/base/` | Base force-field type definitions and SPONGE writer registration for mass, charge, LJ, bond, angle, dihedral, nb14, cmap, virtual atom, UB, RB, SW, EDIP, listed forces, soft bonds. | High |
| `forcefield/amber/` | Amber templates, parameter data, water models, ions, GAFF/GAFF2, protein/nucleic acid/lipid/glycam import side effects. | High |
| `load.py` | File and parameter loaders: PDB, MOL2, frcmod, parmdat, coordinate, rst7, GROMACS ffitp/molitp, PSF, GRO. | High |
| `build.py` | Topology construction and file export: SPONGE input, PDB, MOL2, GRO, MindSPONGE energy adapter. | High except MindSPONGE |
| `process.py` | Modeling operations: geometry constraints, solvation, ion/solvent replacement, peptide builder, optimize, Region/Lattice/crystal building. | High |
| `forcefield/special/` | FEP, GB, minimization/fake parameters and special writer behavior. | Medium |
| `forcefield/charmm/`, `opls/`, `martini/`, `sw/`, `edip/`, `sybyl/` | Non-Amber ecosystems and supporting parameter import scripts. | Medium |
| `__main__.py`, `tools/` | CLI utilities: converter, mask/exclude generation, trajectory analysis, json2csv, name mapping, mol2rfe, MM/GBSA, unit test runner. | Low for modeling API |
| `analysis/` | Post-processing and analysis helpers: MD analysis, SASA, WHAM. | Low |
| `mdrun/` | Wrapper for running external SPONGE executables. | Low / optional |
| `tools/unittests/` | Original test fixtures and regression intent. | Baseline source, not runtime |

## 2. XpongeCPP Current Architecture

XpongeCPP replaces the original nested Python object graph with flat C++ storage
and a thin Python compatibility layer.

```text
XpongeCPP
├── cpp/core/ and cpp/core.hpp
│   Flat data model: Atom, Residue, ResidueType, Molecule, Topology, IDs.
├── cpp/io/
│   PDB, MOL2, PSF, GRO, coordinate, rst7, SPONGE writers.
├── cpp/forcefield/
│   Amber parser/registry, FEP helpers, non-Amber fixture-level parsers.
├── cpp/topology/
│   Bond graph, angle/dihedral/exclude/nb14 topology construction.
├── cpp/assign/
│   C++ Assign core, GAFF/GAFF2 typing, bond order, ring markers, TPACM4.
├── cpp/solvation/
│   Solvent box and ion replacement.
├── cpp/python/
│   pybind11 API surface.
└── src/XpongeCPP/
    Python compatibility layer, forcefield import modules, GROMACS iterator,
    optional RESP/RDKit/PubChemPy helpers, packaged data.
```

### Core data model

| XpongeCPP type | Purpose | Design notes |
| --- | --- | --- |
| `AtomId`, `ResidueId` | Strong integer IDs for cross-array references. | Avoids raw pointers and supports bounds validation. |
| `Atom` | Atom name/type/element, residue ID, PDB metadata, coordinates, charge/mass, special state fields. | Stored in `Molecule.atoms` as a contiguous vector. |
| `Residue` | Residue name/type metadata, PDB chain/resseq metadata, atom range. | Stores `atom_begin` and `atom_count`; no nested atom ownership. |
| `ResidueType` | Writable residue template with atoms, connectivity, head/tail/connect metadata and version counter. | Used by forcefield import and molecule construction. |
| `ResidueLink` | Explicit cross-residue or explicit mol2/PSF/GROMACS bond. | `residue_links` and `explicit_bonds` are separate from atom storage. |
| `Molecule` | Main structure owner: flat atoms, residue ranges, explicit bonds, residue links, box, special force records, optional topology override. | Python `Molecule` views expose residue/atom-like access without restoring old nested ownership. |
| `Topology` | Built bond/angle/dihedral/exclusion/nb14 arrays. | Kept separate from atom/residue objects. |
| `Assign` | Small-molecule assignment graph: elements, coordinates, charges, formal charges, bonds, markers, atom types, rings. | Hot paths for GAFF, bond order, ring markers and TPACM4 are C++. |

### Main data flow

1. Python import layer loads forcefield modules and calls C++ registries.
2. Structure loaders produce a flat `Molecule` with residue ranges and explicit
   connectivity.
3. Modeling operations mutate the flat `Molecule` in batches.
4. `Save_SPONGE_Input` builds or reuses topology arrays and writes files by
   linear scans.
5. Python only wraps API compatibility, optional dependencies and low-frequency
   user customization paths.

## 3. Supported and Missing Features

### Public Python API and namespace

| Feature | Status | Notes |
| --- | --- | --- |
| `import XpongeCPP as Xponge` style | Supported | Common load/build/modeling aliases are exposed. |
| Snake/camel/underscore aliases for common APIs | Partial | Common names are present; original global alias generation is not fully replicated. |
| Original dynamic global namespace injection for all `AtomType`/`ResidueType` names | Partial | Templates are registered and queryable, but Xponge's full global variable behavior is not restored. |
| `GlobalSetting` | Partial | GROMACS include path and bonded parser hooks exist for GROMACS compatibility; full original settings system is not complete. |
| Old Python object model compatibility | Not supported | XpongeCPP intentionally uses C++ flat data structures. |

### Assign

| Feature | Status | Notes |
| --- | --- | --- |
| `Assign` graph object and basic IO entry points | Supported | MOL2/PDB/XYZ/residuetype/SMILES/PubChem/CIF entry points exist. |
| GAFF/GAFF2 assignment from MOL2 | Supported | Rule names and 100-molecule baseline are covered. |
| Bond order, ring markers, kekulize, formal charge | Supported | Main search is C++; custom `penalty_scores` and `extra_criteria` are supported. |
| `AssignRule` custom Python registry | Supported | User custom rules run in Python and do not affect built-in C++ hot paths. |
| TPACM4 | Supported | C++ implementation with regression tests. |
| Gasteiger | Supported with optional dependency | Uses RDKit when installed. |
| RESP | Partial | PySCF backend exists; small-system tests pass. Large and all-parameter Xponge RESP parity is not complete. |
| pH model | Partial | Common phenol/carboxyl/alcohol behavior exists; full original edge coverage is not complete. |
| PubChem real network behavior | Partial | Signature and dependency behavior are present; live network regression is opt-in. |
| CIF symmetry/crystal expansion | Partial | Basic cell/fractional coordinate support exists; full original CIF crystallographic behavior is not complete. |

### Loaders

| Feature | Status | Notes |
| --- | --- | --- |
| `load_pdb` | Supported for Amber/common paths | Chain-local termini, HIS, ACE/NME, OXT, TER, CRYST1, altloc, insertion code, SSBOND/LINK/CONECT, hybrid-36 are covered by tests. |
| `load_mol2` | Supported | Multiple residues, atom type/name/charge/coord, explicit bonds, typed small molecules and custom solvent paths. |
| `load_frcmod` / `load_parmdat` | Supported for Amber | Includes Amber bonded/LJ/CMAP/NB14 paths needed by current workflows. |
| `load_coordinate` | Supported | SPONGE coordinate update and 3/6-value box support. |
| `load_rst7` | Supported | Coordinates and optional box; variants with velocity are not exhaustively covered. |
| `load_gro` | Supported for tested cases | Orthogonal/triclinic box, fixed coordinate columns, velocity ignore, `read_box_angle=False`. |
| `load_molpsf` | Partial | NATOM/NBOND, header variants, connectivity split, charge/type conflict handling. Full PSF ecosystem sections are intentionally skipped like current target behavior. |
| `GromacsTopologyIterator` | Supported | Include stack, macro define/undef, ifdef/ifndef/else/endif, continuation lines, comments. |
| `load_ffitp` | Partial to supported common Xponge behavior | Returns Xponge-style parameter buffers for common atomtypes, pairtypes, bondtypes, angletypes, dihedraltypes, cmaptypes, nonbond params. Unsupported GROMACS function types still raise. |
| `load_molitp` | Partial | Builds molecules/system from moleculetype/atoms/bonds/system/molecules with head/tail prefix. More complex GROMACS molecule blocks and custom bonded parser effects are not fully migrated. |
| Full CHARMM/GROMACS/OPLS topology ecosystem | Partial | Fixture-level non-Amber parsers exist; not yet full Xponge equivalence. |

### Build, topology and export

| Feature | Status | Notes |
| --- | --- | --- |
| Amber topology construction | Supported for current regression set | Uses explicit bonds/templates/residue links and Xponge-like matching semantics for current Amber/GAFF workflows. |
| Core SPONGE input files | Supported | `residue`, `resname`, `atom_name`, `atom_type_name`, `coordinate`, `mass`, `charge`, `LJ`, `bond`, `angle`, `dihedral`, `exclude`, `nb14`. |
| Extra/special SPONGE files | Partial to supported by type | `virtual_atom`, `improper_dihedral`, `cmap`, `nb14_extra`, `urey_bradley`, `Ryckaert_Bellemans`, `listed_forces`, `gb`, `fake_*`, `subsys_division`, `LJ_soft_core`, `SW`, `EDIP`, `bond_soft` have writer paths. Full original special module behavior is not complete. |
| `save_pdb` | Supported for tested Xponge semantics | Chain/TER, CRYST1, SEQRES, SSBOND/LINK/CONECT reconstruction and hybrid-36 have tests. |
| `save_mol2` | Supported for common molecule export | Not a byte-for-byte original Xponge clone. |
| `save_gro` | Supported for tested cases | Box and formatting behavior covered by current fixtures. |
| `build_bonded_force` as old Python API | Not supported | C++ topology builder replaces it internally. |
| Dynamic `Molecule.Set_Save_SPONGE_Input` Python plugin writers | Not supported | Standard writers are built in; arbitrary user writer registration is not replicated. |
| `get_mindsponge_system_energy` | Not planned for core | MindSpore/MindSPONGE runtime adapter is outside current XpongeCPP scope. |

### Process and modeling operations

| Feature | Status | Notes |
| --- | --- | --- |
| `Add_Solvent_Box` | Supported for current workflows | Xponge-style distance scalar/3/6 inputs, arbitrary solvent molecule copying, current TIP3P/Amber regression counts. |
| `Add_Ions` | Supported | Random deterministic replacement by seed, custom solvent residue name. |
| `Set_Box_Padding` | Supported | Explicit box handling and implicit coordinate export behavior are tested. |
| Molecule merge, `+`, `|`, `*` semantics | Supported for tested paths | Link/no-link/repeat behavior covered. |
| `impose_bond`, `impose_angle`, `impose_dihedral` | Not supported | Geometry manipulation routines are not migrated. |
| `h_mass_repartition` | Not supported | No equivalent batch mass repartitioning API yet. |
| `solvent_replace` | Not supported | Ion replacement exists; general selection-based solvent replacement does not. |
| `main_axis_rotate` | Not supported | No principal-axis rotation helper yet. |
| `get_peptide_from_sequence` | Not supported | No sequence-to-peptide builder yet. |
| `sort_atoms_by` | Not supported | No template-based atom reorder helper yet. |
| `optimize` | Not supported | Original external SPONGE minimization helper is not migrated. |
| `Region`, `UnionRegion`, `IntersectRegion`, `BlockRegion`, `SphereRegion`, `FrustumRegion`, `PrismRegion` | Not supported | Region boolean geometry is still missing. |
| `Lattice` and crystal building | Not supported | SC/BCC/FCC/HCP/DIAMOND/custom lattice and periodic bonds are not migrated. |

### Forcefields

| Feature | Status | Notes |
| --- | --- | --- |
| Amber ff14SB | Supported | 1KV2 workflows and packaged data are covered. |
| Amber ff19SB / RSFF2C CMAP | Supported for current tests | CMAP generation and export are tested. |
| Amber nucleic acid/lipid/glycam modules | Partial | Import modules and packaged data exist; broad molecular regression is not complete. |
| TIP3P/TIP4P/TIP4PEW/OPC/SPCE | Supported for current tests | Multi-site virtual atom output covered. |
| GAFF/GAFF2 | Supported for current typed MOL2 and assignment tests | Larger ChEMBL full-run remains optional/local. |
| Common ions | Supported for Amber water models currently tested | Broad ion/water-model matrix still needs larger regression. |
| CHARMM/OPLS/GROMACS parsers | Partial | Minimal/parser-fixture functionality; full ecosystem equivalence is not complete. |
| Martini | Not supported | Data may be packaged in reference tree, but runtime model construction is not migrated. |
| SW/EDIP | Partial | Parameter parser/writer fixture paths exist. |
| SYBYL | Partial | Assign/save typing support exists; full original helper ecosystem not complete. |

### CLI, tools, analysis and run wrappers

| Feature | Status | Notes |
| --- | --- | --- |
| `python -m XpongeCPP` CLI | Not supported | No `XpongeCPP.__main__` equivalent yet. |
| Xponge CLI tools: converter, maskgen, exgen, traj_analysis, json2csv, name2name, mol2rfe, mm_gbsa | Not supported | These are outside current modeling API migration. |
| `analysis/` SASA/WHAM/MD analysis | Not supported | Can be migrated later as Python utilities. |
| `mdrun/` SPONGE executable wrapper | Not supported | Current XpongeCPP focuses on input generation and topology/export. |

## 4. Current Regression Coverage

The current test suite covers the following migration surfaces:

- 1KV2 + ff14SB + TIP3P + NA/CL solvation/export counts and key SPONGE headers.
- B96 typed MOL2 + GAFF/frcmod single molecule and 1KV2+B96 assembly headers.
- PDB chain, terminal, SSBOND/LINK/CONECT, hybrid-36 and round-trip behavior.
- MOL2 explicit bonds and custom solvent/ion replacement.
- Amber data packaging and forcefield import modules, including multi-site waters.
- Assign TPACM4, GAFF rules, optional RDKit/PySCF paths, pH model helpers, custom `AssignRule`, custom bond order parameters.
- PSF/GRO/coordinate/rst7 IO fixtures.
- Non-Amber fixture-level parsers for GROMACS/OPLS/CHARMM/SW/EDIP/FEP softcore.

Typical verification commands:

```bash
rtk uv run pytest -q
pixi run test
pixi run test-assign-full
```

At the time this document was written, the full `rtk uv` suite had `112 passed,
2 skipped`.

## 5. Migration Priorities From Here

Recommended order if the goal is to replace Xponge's modeling API:

1. Finish GROMACS `load_ffitp/load_molitp` edge cases against original Xponge baselines.
2. Migrate `process.py` geometry, Region/Lattice and crystal building.
3. Deepen non-Amber forcefield parsers from fixture-level to Xponge-equivalent.
4. Expand Amber/nucleic acid/lipid/glycam regression beyond import smoke tests.
5. Decide whether dynamic Python writer plugins are required. If not required,
   keep C++ fixed writers as the supported architecture.
6. Treat CLI/tools/analysis/mdrun as a separate compatibility layer after the
   modeling API is stable.
