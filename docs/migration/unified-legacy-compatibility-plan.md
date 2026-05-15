# XpongeCPP Unified Legacy Compatibility Migration Plan

## Goal

Build a single, centralized compatibility layer that lets `XpongeCPP` support the
legacy `Xponge` syntax, import paths, naming style, and common workflow helpers as
fully as practical for core modeling / force-field / IO / process workflows.

This plan replaces the current pattern of scattered one-off shims with a unified
compatibility architecture and an explicit execution process.

## Execution Rule

Before **every** edit round for this migration, re-read this file first:

- [unified-legacy-compatibility-plan.md](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/docs/migration/unified-legacy-compatibility-plan.md)

This file is the source of truth for scope, sequencing, and acceptance.

## Audit Scope

This migration targets the legacy `Xponge` surface used by real scripts and tests in:

- `Xponge/__init__.py`
- `Xponge/helper/__init__.py`
- `Xponge/load.py`
- `Xponge/build.py`
- `Xponge/process.py`
- `Xponge/forcefield/*`
- `Xponge/forcefield/special/*`
- `Xponge/tools/unittests/*`

Out of scope for the first compatibility wave:

- `analysis`
- `cli`
- `mdrun`
- external bundled tools under `tools/`

These may be documented, but are not priority implementation targets for the first
 unified compatibility layer.

To reduce import-time breakage for high-frequency old scripts, the first wave may
still ship import-compatible placeholder shims for these out-of-scope areas so
legacy code fails later with a clear `NotImplementedError` instead of failing
during import.

## What Legacy Xponge Exposes

### 1. Top-level package identity

Legacy users write:

```python
import Xponge
```

not:

```python
import XpongeCPP
```

Many scripts also expect subpackages like:

```python
import Xponge.forcefield.amber.ff19sb
from Xponge.forcefield.special import gb
from Xponge.forcefield.base import bond_base
from Xponge.helper import GlobalSetting
```

### 2. Name-style polymorphism

Legacy `Xponge` intentionally supports multiple names for the same function:

- `load_pdb`
- `Load_Pdb`
- `Load_PDB`
- `LoadPdb`
- `LoadPDB`
- `loadPdb`
- `loadPDB`

This pattern applies broadly to:

- top-level functions
- class methods
- module helpers
- workflow helpers such as `Save_SPONGE_Input`

### 3. Force-field imports have global side effects

Importing a force-field module is not passive. For example:

```python
import Xponge.forcefield.amber.ff14sb
```

is expected to:

- register residue templates
- register force-field parameters
- expose residue symbols like `ALA`, `NALA`, `CALA`
- configure LJ / nb14 / exclude / bonded-type naming behavior

### 4. Main-namespace template injection

After force-field import, users expect:

```python
print(ALA)
mol = NALA + ALA * 10 + CALA
```

without calling `get_template_molecule(...)` manually.

### 5. Legacy algebra / sequence-building syntax

Common legacy expressions include:

```python
NALA + ALA * 10 + CALA
ACE + ALA * 5 + NME
mol = ALA + GLY + SER
```

This implies support for:

- template-like `+`
- template-like `*`
- chained sequence-building expressions
- sequence products that still work with IO and process helpers

### 6. Legacy workflow helper modules

Examples:

```python
from Xponge.forcefield.special import gb
gb.set_gb_radius(mol)
```

Legacy callers expect these helpers at module level, even when the core runtime can
already do the underlying operation elsewhere.

### 7. Legacy object-shape conveniences

Examples already observed in real migration cases:

- `ResidueType.get_type(...).deepcopy(...)`
- `omit_atoms(...)`
- `Xponge.Residue(...)`
- `mol.residue_links = []`
- `mol.residue_links = saved_links`
- `mol.add_residue_link(...)`
- instance-style save methods like `mol.save_pdb(...)`
- top-level helper base names such as `Type`, `Entity`, `AbstractMolecule`, `ResidueLink`
- top-level constants such as `pi`, `kb`, `bar`

## Current XpongeCPP Compatibility Assets

These already exist and must be folded into the new unified layer instead of left as
isolated patches.

### Status Layers

For later completion audit, the current compatibility surface must be interpreted
in three explicit buckets rather than treated as uniformly “done”.

#### 1. Real runnable compatibility

These areas are not only importable or name-compatible; they have direct
repository regression evidence that the legacy call shape actually executes.

Examples currently in this bucket include:

- `Xponge` package-name shim with high-frequency import paths
- top-level IO aliases such as:
  - `Load_PDB / LoadPdb`
- `LoadMol2 / Load_Mol2`
- `SaveMol2`
- `save_gro`
- `Save_GRO`
- `Add_Ions / AddIons`
- `Add_Molecule / AddMolecule`
- `AddSolventBox / Add_Solvent_Box`
- `Impose_Bond / Impose_Angle / Impose_Dihedral`
- `H_Mass_Repartition`
- `Get_Peptide_From_Sequence`
- `SetBoxPadding / Set_Box_Padding`
- `Solvent_Replace`
- `Main_Axis_Rotate`
- `ff19sb + gb` workflow front half, including:
  - `Xponge.NALA / Xponge.ALA / Xponge.CALA`
  - bare-name `NALA + ALA * n + CALA`
  - `gb.set_gb_radius(...)`
  - `gb.SetGbRadius(...)`
  - `Save_PDB(...)`
  - `Save_SPONGE_Input(...)`
- manual 8RYK `frcmod` compatibility
- OPLS/non-amber `poly` compatibility
- `ResidueType` algebra (`+`, `*`)
- `residue_links` API and assignment compatibility
- `Assign` object-shape compatibility for:
  - `atom_types = [...]`
  - `set_atom_type(…, AtomTypeObj)`
  - `to_residuetype(...)` with duplicate-name fallback
  - `deleteBond(...)` / `Delete_Bond(...)`
- first-wave executable `forcefield.base` compatibility for:
  - `LJType.New_From_String(...)` / `get_type(...)`
  - `BondType.New_From_String(...)`
  - `bond_base._gmx_parser(...)`
  - `AngleType.New_From_String(...)`
  - `ProperType.New_From_String(...)`
  - `ImproperType.New_From_String(...)`
  - `ImproperType.same_force(...)`
  - `NB14Type.New_From_String(...)`
  - `VirtualType2.New_From_String(...)`
  - `CMapType.New_From_Dict(...)`
  - `AtomType.New_From_String(...)` interacting with imported
    `mass_base / charge_base / lj_base`
- `save_sponge_input(ResidueType, prefix)` returning a built `Molecule`
- `build_bonded_force(...)` first-wave runnable compatibility for
  `Molecule / Residue / ResidueType`
- real first-wave compatibility for:
  - `analysis.MdoutReader`
  - `analysis.wham.WHAM`
  - `mdrun.run(...)`
- `special.min` fake-parameter export toggles:
  - `save_min_bonded_parameters(...)`
  - `do_not_save_min_bonded_parameters(...)`
  - `Save_Min_Bonded_Parameters(...)`
  - `Do_Not_Save_Min_Bonded_Parameters(...)`
- metadynamics-style front-half workflow
- FEP-style uncovalent and covalent front-half workflows
- executable `special.fep` helper mutation path:
  - `Set_LJ_Type_B(...)`
  - `Set_Subsys(...)`
  - `Enable_LJ_Soft_Core(...)`
  - `Save_Soft_Core_LJ(...)`
  - followed by `Save_SPONGE_Input(...)`
- `Save_GRO(...)` and `mol.Save_GRO(...)` execution
- execution-level `CVSystem` CamelCase aliases

#### 2. Alias / import compatibility only

These areas have evidence that names and import paths exist, but the current
repository evidence is still more shape-level than end-to-end workflow-level.

Examples include:

- broad top-level alias families that are present but not all individually
  exercised in realistic scripts
- many helper namespace convenience exports
- `forcefield.base.*` surfaces beyond the currently exercised first-wave
  registries
- `forcefield.special.fep` helper names beyond the currently exercised workflow
  fronts
- `analysis.md_analysis` and `analysis.sasa` package surfaces when the optional
  `MDAnalysis` dependency is not installed; these now behave as dependency-aware
  conditional shims instead of unconditional placeholders
- `get_mindsponge_system_energy(...)` in environments where the optional
  `mindsponge / mindspore` dependency stack is not installed; the current
  first-wave layer now provides real legacy todo registration and a clear
  dependency-aware `ModuleNotFoundError`, but cannot execute without that stack

#### 3. Import-compatible placeholder shims

These areas are intentionally **not** first-wave runtime-complete, but the
current core-scope migration no longer has a known placeholder-only helper in
its primary completion path.

Audit rule:

- Bucket 1 counts as implementation evidence.
- Bucket 2 counts only as partial evidence.
- Bucket 3 does not count as completed implementation evidence.

### Evidence Matrix

The following repository tests currently provide the strongest direct evidence
for first-wave legacy compatibility:

- [test_compat_layer.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_compat_layer.py)
  covers:
  - package-name shim `Xponge`
  - high-frequency legacy import paths
  - top-level alias call shapes
  - executable `Add_Ions(...)` process workflow alias, including template-like ion keys
  - executable `Add_Molecule(...)` process workflow alias
  - executable `AddSolventBox(...)` process workflow alias
  - executable `Impose_*` and `H_Mass_Repartition(...)` process aliases
  - executable `Get_Peptide_From_Sequence / SetBoxPadding / Solvent_Replace`
    process workflow aliases
  - executable `Main_Axis_Rotate(...)` process alias
  - executable `special.min` fake-parameter export toggles
  - `ResidueType` algebra
  - `ff19sb + gb`
  - `save_mol2 / save_sponge_input / save_gro`
  - `build_bonded_force(...)` first-wave runnable no-op compatibility
  - `Assign` object-shape compatibility
  - metadynamics-style front half
  - FEP-style front halves
  - executable `special.fep` helper mutation path
  - `CVSystem` CamelCase execution aliases
  - first-wave real `MdoutReader` parsing
  - first-wave real `analysis.wham.WHAM`
  - first-wave real `mdrun.run(...)` command-surface behavior

- [test_legacy_import_matrix.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_legacy_import_matrix.py)
  covers:
  - core-scope old `Xponge.*` import paths, including
    `forcefield.special.gb/fep/min`
  - explicit `Xponge.forcefield.base.*` package-name shims for:
    `angle_base`, `charge_base`, `lj_base`, `mass_base`,
    `cmap_base`, `dihedral_base`, `nb14_base`, and `virtual_atom_base`

- [test_8ryk_manual_frcmod.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_8ryk_manual_frcmod.py)
  covers:
  - manual 8RYK `frcmod` compatibility workflow

- [test_poly_non_amber_workflow.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_poly_non_amber_workflow.py)
  covers:
  - non-amber / OPLS `poly` compatibility workflow

- [test_b96_mol2_gaff.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_b96_mol2_gaff.py)
  covers:
  - `save_mol2` connectivity preservation
  - duplicate atom-name MOL2 edge cases
  - B96 / 1KV2 export headers and SPONGE output parity
  - these export checks now execute the `XpongeCPP` workflow in a fresh
    subprocess so they remain stable even when earlier compatibility tests have
    intentionally mutated global legacy force-field state in the parent test
    process

- `tests/test_compat_layer.py::test_forcefield_base_first_wave_minimal_ffbase_compatibility`
  covers:
  - executable first-wave `forcefield.base` compatibility for `LJType`,
    `BondType._gmx_parser`, `AngleType`, `ProperType`, `ImproperType`,
    `ImproperType.same_force(...)`, `NB14Type`, `VirtualType2`, `CMapType`, and
    `AtomType.New_From_String(...)`

Current confirmed command-level evidence from recent runs includes:

- `rtk pixi run pytest -q tests/test_compat_layer.py tests/test_legacy_import_matrix.py`
- `rtk pixi run pytest -q tests/test_8ryk_manual_frcmod.py`
- `rtk pixi run pytest -q tests/test_poly_non_amber_workflow.py`
- `rtk pixi run pytest -q tests/test_b96_mol2_gaff.py`
- `rtk pixi run pytest -q tests/test_compat_layer.py -k ffbase`
- `rtk pixi run pytest -q tests/test_compat_layer.py tests/test_legacy_import_matrix.py tests/test_8ryk_manual_frcmod.py tests/test_poly_non_amber_workflow.py tests/test_b96_mol2_gaff.py`

### Completion Audit Snapshot

Restated current concrete deliverables for this migration:

1. A unified migration / audit document exists and is kept current.
2. A centralized compatibility layer exists under `src/XpongeCPP/_compat`.
3. Old `Xponge` package-name imports and high-frequency submodule imports work.
4. Core first-wave legacy syntax, workflow helpers, and object-shape behaviors
   execute with repository evidence.
5. Remaining unsupported areas are explicitly identified and not miscounted as
   complete.

Current evidence against those deliverables:

- Deliverable 1:
  - this file exists and is updated every implementation round
- Deliverable 2:
  - `_compat` now centralizes alias registration, runtime patch installation,
    assign compatibility, workflow helpers, process-surface wrappers, and
    template symbol syncing
- Deliverable 3:
  - `tests/test_legacy_import_matrix.py`
  - explicit package-name shims under `src/Xponge/`
- Deliverable 4:
  - `tests/test_compat_layer.py`
  - `tests/test_8ryk_manual_frcmod.py`
  - `tests/test_poly_non_amber_workflow.py`
  - `tests/test_b96_mol2_gaff.py`
  - combined green evidence command with 69 passing tests
- Deliverable 5:
  - unsupported / conditional areas are separated into status buckets above

Current first-wave completion conclusion:

- The core modeling / force-field / IO / process compatibility scope defined at
  the top of this document now has repository evidence for:
  - centralized `_compat` installation and alias registration
  - executable `Xponge` package-name shim and high-frequency import paths
  - executable legacy workflow helpers and object-shape compatibility
  - first-wave executable `forcefield.base` registries
  - dependency-aware `get_mindsponge_system_energy(...)` compatibility surface
- Remaining caveats are now dependency-conditional or explicitly outside the
  first implementation wave, rather than unresolved core-scope placeholders.

### A. Top-level aliases and helpers

Current code already exposes many legacy-friendly names from:

- [__init__.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/__init__.py)
- [process.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/process.py)

Including examples such as:

- `Save_PDB`
- `Save_Mol2`
- `Save_GRO`
- `Save_SPONGE_Input`
- `Load_Parameter_From_FFITP`

The current top-level `XpongeCPP` export surface also includes a first-wave set of
legacy helper base names and constants expected by old `Xponge` scripts, including:

- `Atom`
- `Residue`
- `Type`
- `Entity`
- `AbstractMolecule`
- `ResidueLink`
- `AtomType`
- `pi`
- `kb`
- `bar`

The core names imported directly by the legacy `Xponge.__init__` surface are now
present in `XpongeCPP` as first-wave compatibility exports, including:

- assignment helpers such as `Assign` and `get_assignment_from_*`
- helper/runtime names such as `set_global_alternative_names`, `source`,
  `generate_new_pairwise_force_type`, and `generate_new_bonded_force_type`
- build/process names such as `save_*`, `impose_*`, `add_solvent_box`,
  the `Region` family, and `Lattice`
- first-wave runnable compatibility for `build_bonded_force`
- dependency-aware first-wave compatibility for `get_mindsponge_system_energy`

The top-level alias surface is no longer only “name present”; it now also has
direct regression coverage for real old-style I/O invocation forms, including:

- `Load_PDB(...)`
- `LoadPdb(...)`
- `SavePDB(..., write_cryst1=False)`
- `LoadMol2(..., ignore_atom_type=True)`
- `Load_Mol2(..., as_template=True)`
- `SaveMol2(...)`
- `save_mol2(obj)` with legacy default output naming when `filename=None`
- `save_gro(mol, path)`
- `import Xponge.forcefield.amber.ff19sb`
  followed by module-level template globals such as `Xponge.NALA`, `Xponge.ALA`,
  and `Xponge.CALA`
- `from Xponge.forcefield.special import gb` together with
  `Xponge.NALA + Xponge.ALA * n + Xponge.CALA`, `gb.set_gb_radius(mol)`,
  `Save_PDB(...)`, and `Save_SPONGE_Input(...)`

For `load_mol2`, the current first-wave compatibility rule is:

- normal MOL2 loading will try to mirror legacy residue-type registration
- but if template registration fails specifically because the MOL2 carries
  duplicate atom names that are still acceptable for ordinary molecule loading,
  the compatibility layer falls back to plain molecule loading instead of
  failing the read

The helper-side first-wave compatibility also now includes a minimal but usable
`AtomType` registry surface for old tests and scripts that rely on:

- `AtomType.New_From_String(...)`
- `AtomType.get_type(...)`
- `AtomType.Get_Type(...)`
- `AtomType.get_all_types(...)`

For old Amber water workflows, the first-wave `tip3p` import side effects now
also expose the most common atom-type names expected by legacy scripts:

- `AtomType.get_type("HW")`
- `AtomType.get_type("OW")`

The object-surface first-wave compatibility also now includes the most common
old hand-built modeling shapes used by the early legacy tests, including:

- `Assign.addAtom(...)`
- `Assign.Add_Atom(...)`
- `Assign.addBond(...)`
- `Assign.Add_Bond(...)`
- `ResidueType(name="...")`
- `Molecule(name="...")`
- `ResidueType.add_atom(..., AtomTypeObj, ...)`
- `ResidueType.atoms` on manually constructed residue types

The assignment-side first-wave compatibility now also includes the common old
`AssignRule` DSL aliases and decorator shapes, including:

- `AssignRule.addRule(...)`
- `AssignRule.Add_Rule(...)`
- `AssignRule.setPreAction(...)`
- `AssignRule.Set_Pre_Action(...)`
- `AssignRule.setPostAction(...)`
- `AssignRule.Set_Post_Action(...)`

For `Assign` itself, the first wave currently supports common old creation and
typing helpers such as:

- `addAtom(...)`
- `Add_Atom(...)`
- `addBond(...)`
- `Add_Bond(...)`
- `determine_atom_type(...)`
- `Determine_Atom_Type(...)`
- `Determine_Ring_And_Bond_Type(...)`
- `Save_As_PDB(...)`
- `Set_Charge(...)`
- `Set_Charges(...)`
- `Set_Coordinate(...)`
- `Set_Formal_Charge(...)`
- `Set_Atom_Type(...)`
- `Add_Bond_Marker(...)`
- `Has_Bond_Marker(...)`
- `deleteBond(...)`
- `Delete_Bond(...)`

The first-wave `Assign` compatibility now also includes several high-frequency
object-shape behaviors used directly by old workflow tests:

- `assign.set_atom_type(index, AtomTypeObj)`
- `assign.Set_Atom_Type(index, AtomTypeObj)`
- `assign.atom_types = [...]`
- `assign.to_residuetype("NAME")` when atom names would otherwise collide
  (the compatibility layer now falls back to the old `add_index_to_name` style,
  producing names such as `O`, `O1`, `H`, `H1`)

The first-wave SPONGE writer compatibility now also follows the old high-level
Python surface more closely:

- `save_sponge_input(ResidueType, prefix)` now works
- `save_sponge_input(...)` and `Save_SPONGE_Input(...)` now return the built
  `Molecule`, matching legacy `Xponge.build.save_sponge_input`
- `add_solvent_box(ResidueType, WAT, distance, ...)` now works and returns the
  built solvated `Molecule`, matching the common old FEP workflow shape
- a metadynamics-style front-half workflow
  (`get_assignment_from_smiles -> atom_types -> to_residuetype -> save_sponge_input
  -> CVSystem.add_cv_dihedral/meta1d/print/output`) now has direct regression
  coverage
- the front halves of both old FEP workflows now have direct regression
  coverage:
  - uncovalent: `SMILES -> gaff typing -> TPACM4 charge -> to_residuetype ->
    add_solvent_box -> save_pdb/save_mol2`
  - covalent: `ACE + ALA + NME -> add_solvent_box -> save_pdb/save_mol2`

For helper/workflow naming style, the first-wave compatibility now also includes
legacy alias forms on high-frequency CV and GB helpers:

- `CVSystem.Add_Center / AddCenter / addCenter`
- `CVSystem.Add_CV_Position / AddCVPosition / addCvPosition`
- `CVSystem.Add_CV_Dihedral / AddCVDihedral / addCvDihedral`
- `CVSystem.Add_CV_RMSD / AddCVRMSD / addCvRMSD`
- `CVSystem.Print / Steer / Restrain / Meta1D / Meta_1D / Output`
- `gb.SetGbRadius / gb.setGbRadius`

These aliases are no longer only “name present”; there is now direct regression
coverage for a real CamelCase workflow shape using:

- `gb.SetGbRadius(...)`
- `CVSystem.Add_Center(...)`
- `CVSystem.Add_CV_Position(...)`
- `CVSystem.Print(...)`
- `CVSystem.Steer(...)`
- `CVSystem.Output(...)`

There is also direct regression coverage for the CamelCase dihedral-bias path:

- `CVSystem.Add_CV_Dihedral(...)`
- `CVSystem.Restrain(...)`
- `CVSystem.Meta1D(...)`
- `CVSystem.Output(...)`

There is also direct regression coverage for the CamelCase RMSD path:

- `CVSystem.Add_CV_RMSD(...)`
- `CVSystem.Print(...)`
- `CVSystem.Output(...)`

For the `Xponge` package-name shim itself, the first-wave compatibility now also
recreates the old “script-global” style more closely:

- importing `Xponge` injects high-frequency top-level names such as
  `Save_PDB` and `Save_SPONGE_Input` into `__main__`
- force-field imports that register templates now also make names such as
  `NALA`, `ALA`, and `CALA` available in `__main__`
- a bare-name script shape like
  `mol = NALA + ALA * n + CALA; gb.Set_GB_Radius(mol); Save_PDB(...); Save_SPONGE_Input(...)`
  now has direct regression coverage

### B. Existing lightweight compat module

- [compat.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/compat.py)

This currently provides:

- instance-style molecule save methods
- namespace injection of template globals
- public re-export of runtime patch installers and residue-link override helpers
- public re-export of helper-side bonded-force generator compatibility entrypoints

This is useful, but too narrow to serve as the full compatibility architecture.

### B2. Centralized internal compatibility package

- [\_compat](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/__init__.py)

This now exists as the internal consolidation point for:

- import/path shims
- runtime patch installers
- alias registries
- residue-type algebra
- namespace helpers
- workflow helpers

The public surface still includes `compat.py`, but new legacy behavior should be
driven through `_compat` first.

Recent consolidation work now includes:

- [runtime.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/runtime.py)
  as the centralized installer for legacy `ResidueType` / `Molecule` / `Residue`
  object-surface patches
- [assign.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/assign.py)
  as the centralized installer for legacy `Assign` object-shape, DSL, charge,
  and `to_residuetype(...)` compatibility patches that were previously patched
  inline in `XpongeCPP.__init__`
- [api.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/api.py)
  and [compat.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/compat.py)
  now expose `install_legacy_assign_patches(...)` alongside the existing runtime
  and alias installers
- [aliases.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/aliases.py)
  as the centralized registry for top-level legacy alias names instead of a large
  inline map in `XpongeCPP.__init__`
- [imports.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/imports.py)
  now also carries the first-wave helper compatibility implementations for
  `Generate_New_Bonded_Force_Type` and `Generate_New_Pairwise_Force_Type`
- [workflows.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/workflows.py)
  now centralizes high-frequency legacy workflow wrappers and shared state for:
  - `gb.set_gb_radius(...)`
  - `build_bonded_force(...)`
  - `get_mindsponge_system_energy(...)` together with first-wave
    `MindSponge todo` registration and dependency-aware error handling
  - `special.min` bonded-parameter export toggles
- [process.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/process.py)
  now centralizes the high-frequency legacy process-surface wrappers that used
  to live directly in `XpongeCPP.process`, including:
  - `Add_Ions`
  - `Add_Molecule`
  - `Add_Solvent_Box`
  - `Set_Box_Padding`
  - `Save_PDB`
  - `Save_Mol2`
  - `Save_GRO`
  - `Save_SPONGE_Input`
  while `XpongeCPP.process` keeps the lower-level geometry / molecule
  transformation implementations

First-wave import-compatible placeholder bridges now also exist for the most
frequently imported out-of-scope modules:

- [analysis/__init__.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/analysis/__init__.py)
- [analysis/md_analysis.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/analysis/md_analysis.py)
- [analysis/sasa.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/analysis/sasa.py)
- [mdrun.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/mdrun.py)

### C. Legacy residue / molecule behavior patches

- [legacy_types.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/legacy_types.py)

Current compatibility here includes:

- `ResidueType.get_type(...)` legacy handle
- `+`
- `*`
- `deepcopy(...)`
- `omit_atoms(...)`
- `Xponge.Residue(...)`-style conversion path
- `add_residue(...)`
- `residue_links` explicit API and property override

The underlying behavior still lives in `legacy_types.py`, but the installation of
these object-surface patches is now routed through:

- [runtime.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/runtime.py)

`process.Save_PDB(...)` now also consumes residue-link override state through
`_compat.runtime` instead of importing `legacy_types` internals directly.

### D. FRCMOD / GAFF / parmchk2 compatibility

- [forcefield/amber/gaff.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/forcefield/amber/gaff.py)

Current compatibility here includes:

- `parmchk2_gaff(...)`
- path or object input
- `XpongeLib` runtime bridging

### E. OPLS / non-amber compatibility

- [forcefield/opls/__init__.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/forcefield/opls/__init__.py)

Current compatibility here includes:

- `load_parameter_from_ffitp(...)`
- OPLS force-field state registration

### E2. Packaged legacy data-module bridges

Some packaged legacy data modules still expect the old relative import layout.
Current compatibility now includes:

- [data/base/__init__.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/src/XpongeCPP/data/base/__init__.py)
  as a bridge from `XpongeCPP.data.base` to `XpongeCPP.forcefield.base`
- first-wave thin shims for:
  - `cmap_base`
  - `dihedral_base`
  - `nb14_base`
  - `virtual_atom_base`

This is sufficient for `import XpongeCPP.data.amber` to succeed again under the
current compatibility surface.

### F. Current regression assets

Important existing tests that should be treated as compatibility assets:

- [test_compat_layer.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_compat_layer.py)
- [test_8ryk_manual_frcmod.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_8ryk_manual_frcmod.py)
- [test_poly_non_amber_workflow.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_poly_non_amber_workflow.py)
- [test_non_amber_parsers.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_non_amber_parsers.py)
- [test_b96_mol2_gaff.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_b96_mol2_gaff.py)
- [test_pdb_chain_terminal_semantics.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_pdb_chain_terminal_semantics.py)
- [test_xpongecpp_api.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_xpongecpp_api.py)
- [test_legacy_import_matrix.py](/media/yuh/BCDC9249DC91FDB8/Software/Xponge/Xponge-CPP/tests/test_legacy_import_matrix.py)

### G. Package-name and import-path compatibility already landed

There is now a real `src/Xponge/` shim package with working coverage for a first
wave of high-frequency legacy imports, including:

- `Xponge.load`
- `Xponge.build`
- `Xponge.process`
- `Xponge.assign`
- `Xponge.helper`
- `Xponge.helper.file`
- `Xponge.helper.gromacs`
- `Xponge.helper.cv`
- `Xponge.helper.math`
- `Xponge.helper.namespace`
- `Xponge.analysis`
- `Xponge.analysis.md_analysis`
- `Xponge.analysis.sasa`
- `Xponge.mdrun`
- `Xponge.tools`
- `Xponge.tools.unittests`
- `Xponge.forcefield.amber.*` first-wave high frequency modules
- `Xponge.forcefield.base.*` first-wave high frequency modules
- `Xponge.forcefield.charmm.*` high frequency modules
- `Xponge.forcefield.opls.oplsaam`
- `Xponge.forcefield.special.gb`
- `Xponge.forcefield.sw.mw`

As of the current migration round, the high-frequency core-scope legacy imports from
`Xponge/tools/unittests` have been smoke-tested successfully in `XpongeCPP`, including:

- `Xponge`
- `Xponge.forcefield.amber.tip3p`
- `Xponge.forcefield.amber.ff14sb`
- `Xponge.forcefield.base.charge_base`
- `Xponge.forcefield.base.mass_base`
- `Xponge.load`
- `Xponge.forcefield.charmm.tip3p_charmm`
- `Xponge.forcefield.amber.gaff`
- `Xponge.forcefield.amber.rsff2c`
- `Xponge.forcefield.charmm.charmm36`
- `Xponge.forcefield.charmm.charmm27`
- `Xponge.forcefield.amber.bsc1`
- `Xponge.forcefield.sw.mw`
- `Xponge.forcefield.opls.oplsaam`
- `Xponge.forcefield.base.lj_base`
- `Xponge.forcefield.amber.tip4pew`
- `Xponge.forcefield.amber.ol3`
- `Xponge.forcefield.amber.ff19sb`
- `Xponge.build`
- `Xponge.helper`
- `Xponge.helper.file`
- `Xponge.helper.gromacs`
- `Xponge.helper.cv`
- `Xponge.helper.math`
- `Xponge.helper.namespace`
- `Xponge.analysis`
- `Xponge.analysis.md_analysis`
- `Xponge.analysis.sasa`
- `Xponge.mdrun`
- `Xponge.assign`
- `Xponge.tools`
- `Xponge.tools.unittests`

## Known Compatibility Gaps

The following are explicitly known gaps from current audit:

### 1. Unified alias-generation gap

Compatibility names currently live in several files and are not generated from one
registry or one naming policy.

### 2. Centralization gap

Legacy compatibility logic is spread across:

- `__init__.py`
- `process.py`
- `legacy_types.py`
- `compat.py`
- force-field subpackages
- `src/Xponge/*` thin package shims

This makes coverage hard to reason about and easy to regress.

Compared with the start of this migration, two high-churn areas are now partially
centralized:

- object runtime patch installation is routed through `_compat/runtime.py`
- `Assign` legacy monkeypatch installation is routed through `_compat/assign.py`
- top-level alias registration is routed through `_compat/aliases.py`

The remaining work is to keep moving behavior out of scattered module-level patch
sites and into `_compat` installers without regressing the existing scripts.

### 3. Heavy helper semantics gap

Some helper modules can no longer be treated as pure import shims. The most visible
case is:

- `Xponge.helper.cv.CVSystem`
- `Xponge.GlobalSetting` / `Xponge.helper.GlobalSetting`
- `Xponge.helper.source`
- `Xponge.helper.Xopen`

A first-wave minimal compatibility implementation now exists for:

- `add_center(...)`
- `add_cv_position(...)`
- `add_cv_dihedral(...)`
- `add_cv_rmsd(...)`
- `print(...)`
- `steer(...)`
- `restrain(...)`
- `meta1d(...)`
- `output(...)`

Also, common helper-object semantics now exist for:

- `GlobalSetting.purpose`
- `GlobalSetting.logger`
- `GlobalSetting.Add_GMX_Include_Path(...)`
- `GlobalSetting.add_pdb_residue_alias_mapping(...)`
- `GlobalSetting.Set_Invisible_Bonded_Forces(...)`
- `GlobalSetting.Add_Unit_Transfer_Function(...)`
- `source(...)`
- `Xopen(...)`

but this is not yet a full replacement for the original MDAnalysis-backed helper.
Notably, broad MDAnalysis selection compatibility and the full metadynamics helper
surface are still incomplete.

### 4. Out-of-scope but still missing old paths

The following high-level areas remain intentionally outside the first-wave
compatibility target and should not be treated as done:

- `analysis`
- `cli`
- `mdrun`
- `XpongeMoleculeReader`

Current status for these paths is only:

- import-compatible package shims
- placeholder objects or functions
- clear `NotImplementedError` on attempted execution

These may still appear in old scripts, but they are not first-wave acceptance
targets for the current migration.

To reduce import-time breakage in these areas, the current first wave also
includes import-compatible placeholder shims for selected legacy modules such as:

- `Xponge.forcefield.special.min`

These shims are expected to import successfully and fail later with a clear
`NotImplementedError` when their unsupported runtime helpers are actually used.

## Target Architecture

Create a centralized compatibility package under:

- `src/XpongeCPP/compat/`

Recommended submodules:

- `compat/api.py`
  Central public compatibility entrypoints and setup
- `compat/aliases.py`
  Canonical function-alias registry and automatic name generation
- `compat/imports.py`
  Legacy module-path shims and package alias behavior
- `compat/symbols.py`
  Template / residue symbol injection into module namespaces
- `compat/residuetype.py`
  Legacy `ResidueType` handle semantics, algebra, deepcopy, omit, factory helpers
- `compat/molecule.py`
  `residue_links`, instance save methods, legacy object-level behavior
- `compat/forcefield.py`
  Side-effect orchestration for amber/opls/charmm/martini/special imports
- `compat/workflows.py`
  High-level shims like `parmchk2_gaff`, `gb`, `min`, `fep`

### Optional package-level shim

To support true legacy imports, add a separate package shim:

- `src/Xponge/`

This package should be as thin as possible and delegate into `XpongeCPP.compat`.

Examples:

- `src/Xponge/__init__.py`
- `src/Xponge/forcefield/...`
- `src/Xponge/helper/...`

The `Xponge` shim package should:

- forward old imports
- expose legacy names
- avoid reimplementing real logic

All real compatibility logic should stay in `XpongeCPP.compat.*`.

Implementation note:

- the current internal consolidation package is `src/XpongeCPP/_compat/`
- `src/XpongeCPP/compat.py` remains the public stable entrypoint
- when `src/Xponge/__init__.py` mirrors public names from `XpongeCPP`, it must
  avoid overwriting shim-backed subpackages such as `forcefield`, `helper`,
  `load`, `build`, `process`, `assign`, `analysis`, `mdrun`, and `tools`

The long-term target can still evolve toward a fuller `compat/` package surface, but
current work should continue consolidating into `_compat` rather than creating new
scattered shims.

## Compatibility Principles

1. **One canonical implementation**
   Legacy helpers should delegate into one authoritative implementation, not duplicate
   logic in multiple files.

2. **Thin shims**
   Package-path shims and module-level shims should be minimal wrappers around the
   centralized compatibility layer.

3. **State management is explicit**
   Global force-field side effects should be coordinated through one compatibility
   force-field manager.

4. **Tests precede expansion**
   Every new legacy syntax feature should land with a regression test or a workflow
   fixture.

5. **Do not regress already-migrated workflows**
   Existing migration work for:

- manual `frcmod`
- `poly`
- `residue_links`
- hybrid-36

must remain covered during consolidation.

## Migration Phases

### Phase 0: Compatibility inventory normalization

Deliverables:

- map all existing compatibility behavior currently in:
  - `__init__.py`
  - `compat.py`
  - `legacy_types.py`
  - `process.py`
  - force-field helpers
- classify each behavior as:
  - keep
  - move
  - redesign
  - deprecate

Acceptance:

- inventory section in this document remains current

### Phase 1: Central compatibility skeleton

Deliverables:

- create `src/XpongeCPP/compat/` package
- move or wrap current `compat.py` behavior into the new structure
- define one initializer such as:
  - `install_legacy_surface(...)`
  - `install_legacy_aliases(...)`
  - `install_legacy_imports(...)`

Acceptance:

- one obvious place exists for new legacy features

### Phase 2: Alias system unification

Deliverables:

- central registry of canonical names
- generated alternate spellings
- consistent top-level and method-level alias installation

Must cover:

- `load_*`
- `save_*`
- `Add_*`
- `Set_*`
- common CamelCase / mixed-case legacy variants

Acceptance:

- alias creation no longer needs to be hand-scattered across multiple modules

### Phase 3: Main-namespace symbol injection unification

Deliverables:

- one centralized template/global symbol injector
- deterministic refresh on force-field import

Must cover:

- `ALA`
- `NALA`
- `CALA`
- water / ion / ligand template globals

Acceptance:

- force-field imports expose legacy symbols consistently

### Phase 4: ResidueType algebra compatibility

Deliverables:

- support `ResidueType + ResidueType`
- support `ResidueType * int`
- support chained sequence-building expressions

Acceptance:

- this style works:

```python
mol = NALA + ALA * 10 + CALA
```

and can be passed to:

- `Set_GB_Radius`
- `Save_PDB`
- `Save_SPONGE_Input`

### Phase 5: Legacy workflow module shims

Deliverables:

- `forcefield.special.gb`
- `forcefield.special.min`
- `forcefield.special.fep`
- selected `forcefield.base.*` shims

Priority:

1. `gb`
2. `base`
3. `min`
4. `fep`

Acceptance:

- old imports work without changing scripts

### Phase 6: Package alias compatibility

Deliverables:

- optional `src/Xponge/` shim package
- legacy `import Xponge` path support

Acceptance:

- old scripts can import `Xponge` against the compatibility package

### Phase 7: Consolidation of already-implemented compatibility

Deliverables:

- fold existing `frcmod`, `parmchk2`, `residue_links`, `deepcopy`, `omit_atoms`,
  instance save methods, and template namespace injection into the centralized
  compatibility package

Acceptance:

- no important legacy behavior depends on scattered one-off monkey patches alone

## Test Matrix

### A. Alias and import tests

New target file:

- `tests/test_legacy_aliases.py`

Cover:

- top-level aliases
- method aliases
- `Xponge` import-path shims
- `forcefield.special` import-path shims

### B. Namespace and symbol tests

New target file:

- `tests/test_legacy_namespace.py`

Cover:

- force-field import injects `ALA`, `NALA`, `CALA`
- symbol refresh remains correct across multiple force-field imports

### C. ResidueType algebra tests

New target file:

- `tests/test_legacy_residuetype_algebra.py`

Cover:

- `NALA + ALA * 10 + CALA`
- `ACE + GLY * 2 + NME`
- chained save and GB workflows

### D. Workflow shim tests

New target file:

- `tests/test_legacy_workflows.py`

Cover:

- `ff19sb + gb`
- manual `frcmod`
- `poly` non-amber
- representative legacy special/base workflows

### E. Regression preservation

Keep existing files green:

- `tests/test_compat_layer.py`
- `tests/test_8ryk_manual_frcmod.py`
- `tests/test_poly_non_amber_workflow.py`
- `tests/test_non_amber_parsers.py`
- `tests/test_b96_mol2_gaff.py`
- `tests/test_pdb_chain_terminal_semantics.py`
- `tests/test_xpongecpp_api.py`

## Immediate Next Steps

1. Create the centralized `compat/` package skeleton
2. Move current `compat.py` behavior under it
3. Introduce a compatibility registry for aliases and symbol injection
4. Implement `ResidueType` algebra support
5. Add `forcefield.special.gb` shim
6. Add a real `ff19sb + gb` regression against original `Xponge`

## Acceptance for This Planning Stage

This planning stage is complete when:

- this migration document exists
- the document inventories both old `Xponge` semantics and current `XpongeCPP`
  compatibility assets
- the document defines a centralized compatibility architecture
- the document defines phased execution and a reread-before-edit rule

## Notes

- Do not continue adding legacy syntax shims ad hoc in random modules unless the new
  compat layer is explicitly delegating through them.
- When there is a choice between “making this single script pass quickly” and “moving
  the behavior into the centralized compatibility layer”, choose the centralized path.
