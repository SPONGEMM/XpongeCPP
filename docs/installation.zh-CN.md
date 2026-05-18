# XpongeCPP 安装说明（中文）

## 安装目标

当前仓库的打包目标是：

```bash
pip install XpongeCPP
```

安装完成后，同时支持：

```python
import XpongeCPP
import Xponge
```

其中：

- `XpongeCPP` 是新包名
- `Xponge` 是旧脚本兼容包名

当前 wheel 会同时打包：

- `src/XpongeCPP`
- `src/Xponge`

所以旧脚本里的 `import Xponge` 会随着 `XpongeCPP` 一起安装。

## 推荐安装方式

### 1. 普通用户安装

```bash
pip install XpongeCPP
```

安装后建议立即做一次 smoke test：

```bash
python -c "import XpongeCPP, Xponge; print(XpongeCPP.__version__)"
python -c "import Xponge.forcefield.amber.ff19sb; from Xponge.forcefield.special import gb; print('ok')"
```

## 自动依赖策略

基础安装会自动拉取这些高频依赖：

- `numpy`
- `PubChemPy`
- `MDAnalysis`
- `rdkit`
- `pyscf`
- `XpongeLib`

其中 `pyscf` 采用平台条件安装：

- Windows：自动跳过
- Linux / macOS：正常安装

这样可以尽量做到：

- Linux/macOS 用户直接安装后即可使用更多完整功能
- Windows 用户不会因为 `pyscf` 不可用而让整个包安装失败
- `gaff.parmchk2_gaff(...)` 这类 legacy workflow 能自动解析 `XpongeLib` bridge

## 开发安装

如果你是在仓库中开发，推荐继续使用 `pixi`：

```bash
pixi run install-dev
pixi run test
```

也可以直接 editable 安装：

```bash
pip install -e .
```

## 安装后建议验证的典型用法

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

## 构建与发布前验证

仓库里提供了一个本地 PyPI 构建脚本：

```bash
pixi run python scripts/build_pypi.py
```

它会自动完成：

1. 构建 sdist 和 wheel
2. 执行 `twine check`
3. 创建临时虚拟环境
4. 安装构建出的 wheel
5. 执行 `import XpongeCPP` 和 `import Xponge` 的 smoke test

## 如果安装失败

优先检查这些问题：

1. Python 版本是否满足 `>=3.10`
2. 本地是否具备编译 C++ 扩展所需工具链
3. 是否是某个化学依赖在当前平台上不可用

如果你的目标是快速开发和验证，而不是立刻走纯 pip 环境，优先推荐 `pixi`。
