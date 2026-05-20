# XpongeCPP 安装说明

## 推荐安装方式

对于普通用户，目标是：

- `pip install XpongeCPP`
- 安装后同时支持：
  - `import XpongeCPP`
  - `import Xponge`

当前 wheel 会同时打包：

- `src/XpongeCPP`
- `src/Xponge`

因此旧脚本中的 `import Xponge` 兼容层会随包一起安装。

## PyPI 安装

```bash
pip install XpongeCPP
```

安装完成后，建议先做 smoke test：

```bash
python -c "import XpongeCPP, Xponge; print(XpongeCPP.__version__)"
python -c "import Xponge.forcefield.amber.ff19sb; from Xponge.forcefield.special import gb; print('ok')"
```

## 自动依赖策略

基础安装会自动安装以下高频依赖：

- `numpy`
- `PubChemPy`
- `MDAnalysis`
- `rdkit`
- `pyscf`
- `XpongeLib`

其中 `pyscf` 采用平台条件安装：

- Windows 下自动跳过
- 非 Windows 平台正常安装

这意味着：

- Linux / macOS 用户默认会得到完整一些的化学后端
- Windows 用户不会因为 `pyscf` 不可用而导致整包安装失败
- `gaff.parmchk2_gaff(...)` 这类 legacy workflow 能自动解析 `XpongeLib` bridge

## RESP 后端说明

当前 RESP 电荷计算支持多后端策略：

- 默认后端：`PySCF`
- 可选后端：`Psi4`

推荐使用方式：

- Linux / macOS：默认直接使用 `PySCF`
- Windows：安装 `Psi4`，并在 RESP 调用时显式指定 `backend="psi4"`

示例：

```python
assign.calculate_charge("resp", backend="pyscf")
assign.calculate_charge("resp", backend="psi4")
```

如果是在 Windows 上做 RESP，推荐安装顺序是：

```bash
conda install -c conda-forge psi4
pip install XpongeCPP
```

说明：

- `XpongeCPP` 仍然可以正常通过 PyPI 分发
- `Psi4` 被当作可选外部后端，而不是 PyPI 强制依赖
- RESP 数值核心会逐步迁移到 C++，但后端选择接口保持稳定

## 开发安装

如果你是在仓库里开发，推荐继续使用 `pixi`：

```bash
pixi run install-dev
pixi run test
```

也可以用 editable 模式：

```bash
pip install -e .
```

## 安装后建议验证的功能

### 新接口

```python
import XpongeCPP
import XpongeCPP.forcefield.amber.ff19sb
from XpongeCPP.forcefield.special import gb
```

### 旧脚本兼容接口

```python
import Xponge
import Xponge.forcefield.amber.ff19sb
from Xponge.forcefield.special import gb

mol = Xponge.NALA + Xponge.ALA * 2 + Xponge.CALA
gb.set_gb_radius(mol)
```

## 如果安装失败

优先检查：

1. Python 版本是否为 `>=3.10`
2. 是否具备本地编译 C++ 扩展所需工具链
3. 是否是某个可选化学依赖在当前平台不可用

如果只是想快速开发验证，也可以先使用仓库提供的 `pixi` 环境。
