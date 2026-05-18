# XpongeCPP 接口概览（中文）

本文档面向当前第一波可用接口，重点说明：

- 新接口如何使用
- 旧 `Xponge` 脚本如何兼容运行
- 哪些模块是高频主入口

## 两条主入口

### 1. 新入口

```python
import XpongeCPP
```

这是推荐给新代码的用法。

### 2. 旧入口兼容

```python
import Xponge
```

这是为旧脚本保留的兼容包名入口。  
安装 `XpongeCPP` 后，wheel 会同时带上这一层 shim。

## 高频顶层对象

当前高频核心对象包括：

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

高频入口包括：

- `load_pdb`
- `load_mol2`
- `Save_PDB`
- `Save_Mol2`
- `Save_GRO`
- `Save_SPONGE_Input`
- `Add_Solvent_Box`
- `Add_Ions`
- `Set_Box_Padding`

## 旧语法兼容重点

### 模板代数语法

```python
mol = NALA + ALA * 10 + CALA
```

### 旧命名风格

同一函数通常支持多种旧写法，例如：

- `load_pdb`
- `Load_PDB`
- `LoadPDB`

### 旧包路径

例如下面这些旧写法仍然是兼容目标：

```python
import Xponge.forcefield.amber.ff19sb
from Xponge.forcefield.special import gb
from Xponge.helper.file import pdb_filter
```

## 典型工作流

### 模板建模主链

```python
import Xponge
import Xponge.forcefield.amber.ff19sb
from Xponge.forcefield.special import gb

mol = Xponge.NALA + Xponge.ALA * 10 + Xponge.CALA
gb.set_gb_radius(mol)
Xponge.Save_PDB(mol, "ALA.pdb")
Xponge.Save_SPONGE_Input(mol, "ALA")
```

### Assign 判型主链

```python
assign = Xponge.get_assignment_from_smiles("OO")
assign.determine_atom_type("gaff")
res = assign.to_residuetype("TES")
```

## 更适合看源码的模块

如果你想继续理解实现层次，可以优先看：

- `src/XpongeCPP/__init__.py`
- `src/XpongeCPP/_compat/*`
- `src/XpongeCPP/load.py`
- `src/XpongeCPP/build.py`
- `src/XpongeCPP/process.py`
- `src/XpongeCPP/assign/__init__.py`

## 相关文档

- [README.md](../README.md)
- [installation.md](./installation.md)
- [installation.zh-CN.md](./installation.zh-CN.md)
- [api-overview.md](./api-overview.md)
- [xponge-vs-xpongecpp-architecture-status.md](./xponge-vs-xpongecpp-architecture-status.md)
