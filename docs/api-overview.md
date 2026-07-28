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
