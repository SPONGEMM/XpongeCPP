# XpongeCPP

XpongeCPP is a new C++ implementation of common Xponge workflows with a thin Python compatibility layer.

The original Xponge repository is used only as a reference implementation and regression baseline. This repository does not share the original Python object model internally.

## v1 Scope

- Amber-first force-field workflows.
- Flat C++ storage with `Molecule`, `Residue`, and `ResidueType` view semantics.
- Common Python entry points such as `load_pdb`, `load_mol2`, `Add_Solvent_Box`, `Set_Box_Padding`, `Save_SPONGE_Input`, and `Assign`.
- Numeric equivalence goals for SPONGE input, not byte-for-byte compatibility.

## Development Environments

The lightweight development path still uses `rtk uv`:

```bash
rtk uv pip install -e . --force-reinstall --no-cache-dir
rtk uv run pytest -q
```

Full Assign validation needs optional chemistry backends. Use `pixi` for a reproducible environment with RDKit, PubChemPy, and PySCF:

```bash
pixi run install-dev
pixi run test-assign-full
pixi run test-resp
pixi run test
```

PubChem network tests remain opt-in through the test environment; default Pixi tests use local mocks or dependency checks and do not require live network access.
