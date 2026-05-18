# XpongeCPP

## Docs / 文档入口

- English:
  - [Installation Guide](./docs/installation.md)
  - [API Overview](./docs/api-overview.md)
- 中文：
  - [安装说明](./docs/installation.zh-CN.md)
  - [接口概览](./docs/api-overview.zh-CN.md)

XpongeCPP is a new C++ implementation of common Xponge workflows with a thin Python compatibility layer.

The original Xponge repository is used only as a reference implementation and regression baseline. This repository does not share the original Python object model internally.

## Installation

### PyPI-style installation target

The packaging target is:

```bash
pip install XpongeCPP
```

After installation, both of these should work:

```python
import XpongeCPP
import Xponge
```

The wheel is configured to include both:

- `src/XpongeCPP`
- `src/Xponge`

so old `Xponge` package-name imports can continue to work after installation.

### Automatic optional chemistry dependencies

The package now declares a practical default dependency set for pip users:

- `numpy`
- `PubChemPy`
- `MDAnalysis`
- `rdkit`
- `pyscf` on non-Windows platforms
- `mokda-xpongelib`

Windows automatically skips `pyscf` through environment markers.

`mokda-xpongelib` is included so legacy `gaff.parmchk2_gaff(...)` workflows can
resolve the `XpongeLib` bridge automatically after installation.

See:
- [docs/installation.md](./docs/installation.md)

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

## Documentation

- Installation guide:
  - [docs/installation.md](./docs/installation.md)
  - [docs/installation.zh-CN.md](./docs/installation.zh-CN.md)
- API overview:
  - [docs/api-overview.md](./docs/api-overview.md)
  - [docs/api-overview.zh-CN.md](./docs/api-overview.zh-CN.md)
- Architecture / migration status:
  - [docs/xponge-vs-xpongecpp-architecture-status.md](./docs/xponge-vs-xpongecpp-architecture-status.md)
