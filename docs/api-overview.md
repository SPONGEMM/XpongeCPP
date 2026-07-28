# XpongeCPP 接口概览

本文档是当前第一波用户接口说明，重点说明：

- 新接口如何使用
- 旧 `Xponge` 兼容接口如何使用
- 哪些模块是主入口

## 两条主入口

### 1. 新入口

```python
import XpongeCPP
```

推荐给新代码使用。

### 2. 旧入口兼容

```python
import Xponge
```

这是为旧脚本保留的兼容包名入口。

安装 `XpongeCPP` 后，wheel 会同时带上这个 shim 包。

## 顶层高频对象

可直接使用的核心对象包括：

- `Atom`
- `ResidueType`
- `Residue`
- `Molecule`
- `Assign`
- `AssignRule`

## 高频模块

### Force field

```python
import Xponge.forcefield.amber.ff14sb
import Xponge.forcefield.amber.ff19sb
import Xponge.forcefield.amber.gaff
import Xponge.forcefield.amber.tip3p
```

### Special workflows

```python
from Xponge.forcefield.special import gb
from Xponge.forcefield.special import fep
```

### Constrained RESP

`Assign.calculate_charge("RESP", ...)` now accepts general linear constraints
`Cq=d` in addition to equivalence groups:

```python
assignment.calculate_charge(
    "RESP",
    charge=0,
    extra_equivalence=[[1, 2]],
    constraint_matrix=[[1.0, 0.0, 0.0]],
    constraint_targets=[-0.4],
    return_diagnostics=True,
)
diagnostics = assignment.charge_fit_diagnostics
```

Constraints are checked before the QM backend runs. Diagnostics include
constraint rank, dependent rows, KKT condition/residual, ESP error metrics,
and the applied row ledger. Both unconstrained fitting and constrained KKT
solves use the C++ numerical core by default; `core="python"` remains
available for cached-ESP parity checks.

### Load / Build / Process

常见高频入口包括：

- `load_pdb`
- `load_mol2`
- `Save_PDB`
- `Save_Mol2`
- `Save_GRO`
- `Save_SPONGE_Input`
- `save_sponge_input(..., format="raw" | "bundle")`
- `save_sponge_input_raw`
- `save_sponge_input_bundle` (native HighFive/HDF5 topology, protocol, and restart output)
- `Add_Solvent_Box`
- `Add_Ions`
- `Set_Box_Padding`

### Bundle 读取与校验

`XpongeCPP.io_bundle` 提供与 Xponge-origin 同名的 bundle case 和 reader
边界：

```python
from XpongeCPP.io_bundle import BundleReader, bundle_case_from_prefix

case = bundle_case_from_prefix("inputs", "system")
with BundleReader(case) as reader:
    atom_count = reader.read_scalar(
        "topology.spgt.h5", "/topology/atom_count"
    )
```

Reader 会校验 v2 schema、三件套 UUID、topology/atom-order/protocol lineage、
原子维度与 restart 完成状态。已有 bundled mdin 的情况可使用
`scan_bundle_case(...)`。

如需直接分析 bundle 中的 topology 与 restart：

```python
from XpongeCPP.analysis import load_bundle_universe

universe = load_bundle_universe("inputs/system_topology.spgt.h5")
```

该入口返回带原子、残基、键、坐标和周期盒的 MDAnalysis `Universe`。当前
版本会将选择的单帧 materialize 到内存；多帧 trajectory streaming 属于后续
对齐阶段。

受支持的普通力场 bundle 也可以安全地回转为 SPONGE direct/legacy 输入：

```python
from XpongeCPP.io_bundle import convert_bundle_to_legacy

manifest = convert_bundle_to_legacy(
    "inputs", "legacy-inputs", prefix="system"
)
```

转换器先规划并校验所有目标，再原子写入文件；`dry_run=True` 只做规划，
`overwrite=False` 默认拒绝覆盖。当前支持普通 LJ、键、角、周期二面角、
harmonic improper、排除表和两列 nb14。soft-core LJ、CMAP、GB、虚拟原子、
Urey-Bradley、soft bond 与 custom listed force 在 strict 模式下会明确报错，
不会被静默丢弃。

已有 direct/legacy SPONGE case 可以转换为 canonical v2 三件套：

```python
from XpongeCPP.io_bundle import convert_legacy_to_bundle

manifest = convert_legacy_to_bundle(
    "legacy-case", "converted-case", mdin="mdin.spg.toml"
)
```

输出位于 `converted-case/bundle/`，包括 topology、protocol、restart 与
`mdin.bundled.spg.toml`；转换时会重新生成 UUID、canonical content hash 和
lineage。当前只接受上面列出的普通力场子集，遇到其他已绑定输入会明确拒绝。

## 旧语法兼容重点

### 模板代数语法

```python
mol = NALA + ALA * 10 + CALA
```

### 旧命名风格

同一个函数通常支持多种命名：

- `load_pdb`
- `Load_PDB`
- `LoadPDB`

### 旧包路径

例如：

```python
import Xponge.forcefield.amber.ff19sb
from Xponge.forcefield.special import gb
from Xponge.helper.file import pdb_filter
```

## 当前更适合看源码的模块

如果想理解实现层次，可以优先看：

- `src/XpongeCPP/__init__.py`
- `src/XpongeCPP/_compat/*`
- `src/XpongeCPP/process.py`
- `src/XpongeCPP/load.py`
- `src/XpongeCPP/build.py`
- `src/XpongeCPP/assign/__init__.py`

## 相关文档

- [README.md](../README.md)
- [installation.md](./installation.md)
- [xponge-vs-xpongecpp-architecture-status.md](./xponge-vs-xpongecpp-architecture-status.md)
