# XpongeCPP 迁移计划：跑通 8RYK `manual_8ryk_test.py`

## 背景

目标是让 `XpongeCPP` 能直接支持如下 legacy 风格脚本：

1. `import XpongeCPP as Xponge`
2. `import XpongeCPP.forcefield.amber.ff14sb`
3. `import XpongeCPP.forcefield.amber.gaff as gaff`
4. `Xponge.load_mol2(..., as_template=True)`
5. 构造 `PHE/ASP/CCS/GLY` 混合分子
6. 添加 3 条 `residue_link`
7. `Xponge.Save_Mol2(...)`
8. `gaff.parmchk2_gaff(...)`
9. 生成 `manual.raw.frcmod`

当前已确认：

- 旧 `Xponge-origin` 在安装正确的 `XpongeLib` wheel 后，这份脚本可以跑通
- `XpongeCPP` 当前还没有完整支持这条脚本所依赖的 legacy API 面
- `XpongeCPP` 的 `load_mol2(..., as_template=True)` 与 `Save_Mol2(...)` 基础能力基本可用
- 真正缺口主要集中在 `parmchk2_gaff` 兼容入口，以及 `ResidueType` 的 legacy 编辑能力

## 当前现状与缺口

### 已具备的能力

- `XpongeCPP.load_mol2(..., as_template=True)` 已支持模板注册
- `XpongeCPP.Save_Mol2(...)` 已存在，可导出混合体系 mol2
- `Molecule.add_residue_link(...)` 已支持 legacy 风格传入 atom 对象
- `Molecule.residue_links` 已可读取
- `Molecule.atom_index[...]` 兼容代理已存在

### 已确认缺失的能力

- `XpongeCPP.forcefield.amber.gaff` 目前没有正式暴露 `parmchk2_gaff(...)`
- `XpongeCPP._core` 不提供 `_parmchk2`
- `XpongeCPP.forcefield.base` 当前不存在，完整 `frcmod_init.txt` 仍缺这层兼容
- `ResidueType.get_type("PHE")` 当前没有 legacy 风格：
  - `deepcopy(...)`
  - `omit_atoms(...)`
- `molecule.residue_links` 当前是只读属性，不能像旧脚本那样整体赋值替换

### 结论

要跑通 `manual_8ryk_test.py`，不能只补 `parmchk2_gaff` 一项。最少还需要补齐：

1. `gaff.parmchk2_gaff(...)` 公共兼容入口
2. `ResidueType.deepcopy(...)`
3. `ResidueType.omit_atoms(...)`

`forcefield.base` 和完整 boundary 元数据流程不作为本轮必需交付项。

## 实施原则

- 本轮只以“跑通 `manual_8ryk_test.py` 风格脚本”为验收目标
- 不把完整 `frcmod_init.txt` 兼容一起做完
- 不重写现有 mol2 导入导出逻辑，除非新测试证明有实际问题
- `parmchk2` 继续依赖外部 `XpongeLib`，本轮不在 `XpongeCPP` 内重写 `_parmchk2`

## 分步执行计划

### Phase 1：补齐 `gaff.parmchk2_gaff(...)` 公共接口

目标：

- 在 `XpongeCPP.forcefield.amber.gaff` 中正式提供 `parmchk2_gaff(ifname, ofname, direct_load=True, keep=True)`

修改要点：

- 将现有 legacy 风格 `parmchk2_gaff` 逻辑迁入或重新导出到 `src/XpongeCPP/forcefield/amber/gaff.py`
- 运行时从已安装的 `XpongeLib` 包中导入 `_parmchk2`
- 支持两种输入：
  - 直接传 mol2 文件路径
  - 传 `Molecule` / `Residue` / `ResidueType` 风格对象时，先临时导出 mol2
- `direct_load=True` 时，自动 `amber.load_parameters_from_frcmod(ofname, prefix=False)`
- `keep=False` 时，加载后删除输出 frcmod

失败模式要求：

- 如果没有安装 `XpongeLib`，抛出清晰 `ImportError`
- 错误信息要明确指出需要 `mokda-xpongelib` / `XpongeLib`

### Phase 2：补齐 `ResidueType` 的 legacy 编辑能力

目标：

- 让下面这类脚本语义在 `XpongeCPP` 中成立：

```python
phe_type = Xponge.ResidueType.get_type("PHE").deepcopy("PHE_1")
phe_type.omit_atoms([...], charge=None)
```

修改要点：

- 为 `ResidueType` 增加 legacy 风格 `deepcopy(name)` 兼容方法
- 为 `ResidueType` 增加 `omit_atoms(atom_names, charge=None)` 兼容方法

行为要求：

- `deepcopy(name)` 返回一个可独立编辑的新模板类型
- `omit_atoms(...)` 至少要正确支持本 8RYK 脚本里的用法：
  - 根据原模板删去不在 `present` 集合内的原子
  - 不破坏保留原子的类型、坐标、连接信息
  - `charge=None` 表示不强制重分配总电荷策略

本阶段验收重点：

- `PHE/ASP/GLY` 这三种标准残基能按脚本方式裁剪
- 后续 `Xponge.Residue(new_type, directly_copy=True)` 能继续工作

### Phase 3：增加手工 8RYK 脚本回归测试

新增一个精确模拟 `manual_8ryk_test.py` 的测试。

测试内容：

- 导入 `XpongeCPP as Xponge`
- 导入 `XpongeCPP.forcefield.amber.ff14sb`
- 导入 `XpongeCPP.forcefield.amber.gaff as gaff`
- `Xponge.load_mol2(CCS_3.gaff.mol2, as_template=True)`
- 构造 `PHE/ASP/CCS/GLY`
- 对 `PHE/ASP/GLY` 做 `deepcopy + omit_atoms`
- 添加 3 条 `residue_link`
- `Xponge.Save_Mol2(...)`
- `gaff.parmchk2_gaff(...)`

断言：

- `manual.raw.mol2` 成功生成
- `manual.raw.frcmod` 成功生成
- frcmod 非空
- 不出现段错误
- 允许当前 `CO` warning 存在，不把 warning 当作失败

环境处理：

- 如果没有安装 `XpongeLib`，测试应 `skip`
- 一旦环境满足，测试必须是正式通过，不使用 `xfail`

### Phase 4：补小粒度兼容测试，锁住 API 行为

新增两个更小的测试：

1. 文件路径输入测试  
   验证 `parmchk2_gaff(path, out, direct_load=False, keep=True)` 可直接工作

2. 对象输入测试  
   验证 `parmchk2_gaff(molecule_obj, out, ...)` 会先导出临时 mol2，再生成 frcmod

另外补充一个 `ResidueType` 测试：

- `deepcopy + omit_atoms` 后，保留原子集合符合预期

### Phase 5：回归验证与范围控制

至少运行以下测试：

1. 新增 manual 8RYK 回归测试
2. `tests/test_xpongecpp_api.py`
3. `tests/test_b96_mol2_gaff.py`
4. `tests/test_8ryk_regression.py`

预期：

- `spg_init` 继续通过
- `frcmod_init` 仍可保持 `xfail`
- 不能引入现有 GAFF / mol2 / residue link 行为回归

## 非本轮范围

以下内容本轮不要求完成：

- `XpongeCPP.forcefield.base` 全量兼容
- `angle_base` / `bond_base` / `dihedral_base` 旧注册系统
- 完整 `frcmod_init.txt` payload/boundary metadata 工作流
- 在 `XpongeCPP` 内部重写 `_parmchk2`
- 去除 `XpongeLib` 外部依赖

## 验收标准

完成后必须满足：

- `XpongeCPP` 版本的 `manual_8ryk_test.py` 可以成功运行
- 能生成 `manual.raw.mol2` 和 `manual.raw.frcmod`
- 不再出现段错误
- `spg_init` 回归继续通过
- 完整 `frcmod_init` 可以暂时继续保留 `xfail`

## 当前支持情况答复

### 是否支持 omit 原子？

当前 **不支持完整 legacy 用法**。

已确认：

- `ResidueType.get_type("PHE")` 目前没有 `deepcopy(...)`
- `ResidueType.get_type("PHE")` 目前没有 `omit_atoms(...)`

因此你这份 8RYK 手工脚本里对 `PHE/ASP/GLY` 的裁剪逻辑，目前在 `XpongeCPP` 里还跑不起来，必须先补。

### 是否支持修改 residue link？

当前 **部分支持**。

已确认支持：

- `molecule.add_residue_link(atom_a, atom_b)` 可用
- 可以传 atom 对象，兼容层会自动转 index
- `molecule.residue_links` 可读取

当前不支持或不稳定的点：

- `molecule.residue_links = [...]` 这种整体替换目前不行
- `molecule.residue_links = []` 这种清空式写法目前也不行，因为它是只读属性

所以如果后续流程需要“删除、重置、整体覆盖 residue_links”，还需要额外设计一个显式 API，例如：

- `clear_residue_links()`
- `remove_residue_link(...)`
- `set_residue_links(...)`

但对 `manual_8ryk_test.py` 这份脚本本身来说，只要求“追加 3 条 residue link”，当前能力已经够用。
