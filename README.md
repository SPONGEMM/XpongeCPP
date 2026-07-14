# XpongeCPP

> Preview beta: XpongeCPP is currently a preview release and is under active beta testing.
>
> 预览版说明：XpongeCPP 当前仍处于预览阶段，正在进行 Beta 测试。

## Docs / 文档入口

- English:
  - [Installation Guide](./docs/installation.md)
  - [API Overview](./docs/api-overview.md)
  - [Release Guide](./docs/releasing.md)
- 中文：
  - [安装说明](./docs/installation.zh-CN.md)
  - [接口概览](./docs/api-overview.zh-CN.md)
  - [发布说明](./docs/releasing.zh-CN.md)

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
- `geometric`
- `matplotlib`
- `PubChemPy`
- `MDAnalysis`
- `rdkit`
- `pyscf` on non-Windows platforms
- `mokda-xpongelib` (imported as `XpongeLib`)

Windows automatically skips `pyscf` through environment markers.

`mokda-xpongelib` is included so legacy `gaff.parmchk2_gaff(...)` workflows can
resolve the `XpongeLib` bridge automatically after installation.

See:
- [docs/installation.md](./docs/installation.md)

### GitHub Actions packaging matrix

The repository CI currently builds packages on:

- Linux x64: `ubuntu-24.04`
- Linux arm64: `ubuntu-24.04-arm`
- macOS Intel: `macos-15-intel`
- macOS arm64: `macos-15`
- Windows x64: `windows-2025`

Validation is split into two layers:

- all platforms run a minimal wheel smoke test with `numpy` installed
- Linux x64 additionally runs a full dependency install smoke test

This keeps wheel validation broad across operating systems and CPU architectures
without making every matrix job depend on the full optional chemistry stack.

## v1 Scope

- Amber-first force-field workflows.
- Flat C++ storage with `Molecule`, `Residue`, and `ResidueType` view semantics.
- Common Python entry points such as `load_pdb`, `load_mol2`, `Add_Solvent_Box`, `Set_Box_Padding`, `Save_SPONGE_Input`, and `Assign`.
- Numeric equivalence goals for SPONGE input, not byte-for-byte compatibility.

### Amber lipid force fields

Both Lipid17 and Lipid21 are packaged. Import exactly one lipid base family:

```python
import XpongeCPP.forcefield.amber.ff14sb
import XpongeCPP.forcefield.amber.gaff2
import XpongeCPP.forcefield.amber.lipid21  # or lipid17
```

Either lipid import automatically registers the shared PI/phosphoinositide/LysoPL
extension. Input remains Amber's split-residue representation, such as
`PA + SPM + SA`; full lipid names such as `PSM` and `POPC` are not automatically
split. The extension's mixed Lipid/GLYCAM/phosphate/GAFF2 provenance is reported
when it is loaded.

Within one Python process, `ff14sb`/`ff19sb`, `gaff`/`gaff2`, and
`lipid17`/`lipid21` are mutually exclusive pairs. Different families may be
combined; use separate processes to compare alternatives in the same family.

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

RESP supports a multi-backend strategy:

- default backend: `PySCF`
- optional backend: `Psi4`
- shared dispatch layer: `XpongeCPP.qm`

Example:

```python
assign.calculate_charge("resp", backend="pyscf")
assign.calculate_charge("resp", backend="psi4")

from XpongeCPP import qm

qm.run_scf(assign, backend="pyscf")
qm.optimize_geometry(assign, backend="pyscf")
```

For Windows RESP workflows, install `Psi4` separately and select it explicitly:

```bash
conda install -c conda-forge psi4
pip install XpongeCPP
```

PubChem network tests remain opt-in through the test environment; default Pixi tests use local mocks or dependency checks and do not require live network access.

## Documentation

- Installation guide:
  - [docs/installation.md](./docs/installation.md)
  - [docs/installation.zh-CN.md](./docs/installation.zh-CN.md)
- API overview:
  - [docs/api-overview.md](./docs/api-overview.md)
  - [docs/api-overview.zh-CN.md](./docs/api-overview.zh-CN.md)
- Release guide:
  - [docs/releasing.md](./docs/releasing.md)
  - [docs/releasing.zh-CN.md](./docs/releasing.zh-CN.md)
- Architecture / migration status:
  - [docs/xponge-vs-xpongecpp-architecture-status.md](./docs/xponge-vs-xpongecpp-architecture-status.md)

## Packaging roadmap note

The current repository uses a hand-written GitHub Actions workflow plus
`scripts/build_pypi.py` for packaging validation. We are deliberately keeping
this simpler workflow for now because it makes the `Xponge`/`XpongeCPP` dual
package layout and the minimal-vs-full smoke split easy to audit.

`cibuildwheel` is still a good future option once the wheel matrix and release
policy stabilize further, but it is not the current default.
