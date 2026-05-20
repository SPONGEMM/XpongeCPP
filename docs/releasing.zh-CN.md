# XpongeCPP 发布说明（中文）

本文档说明当前仓库如何发布到 PyPI。

## 当前发布策略

仓库现在使用两套 GitHub Actions workflow：

- `build-packages.yml`
  - 用于 push / pull request 验证
  - 在支持的平台矩阵上构建 wheel
  - 在 Linux x64 上构建 sdist
  - 执行元数据检查与 smoke test
- `publish-pypi.yml`
  - 在 GitHub Release 发布后触发
  - 重新构建所有 wheel 和 sdist
  - 通过 Trusted Publishing 上传到 PyPI

## 当前平台矩阵

发布 wheel 当前覆盖：

- Linux x64
- Linux arm64
- macOS Intel
- macOS arm64
- Windows x64

同时在 Linux x64 上构建 sdist。

## Smoke 测试分层

每个 wheel job 都运行一条最小 smoke test：

1. 创建全新虚拟环境
2. 安装 `numpy`
3. 使用 `--no-deps` 安装刚构建出的 wheel
4. 验证：
   - `import XpongeCPP`
   - `import Xponge`

这样可以在不要求所有 runner 都装齐完整化学依赖的前提下，保持 wheel
矩阵的跨平台覆盖。

更完整的依赖安装验证，仍然保留在常规 packaging 验证里，并由 Linux x64
主平台承担。

## PyPI Trusted Publishing 配置

当前发布 workflow 设计为使用 GitHub Actions + PyPI Trusted Publishing。

你需要在 PyPI 项目页里配置 Trusted Publisher，参数应为：

- owner: `yuhaosimba`
- repository: `XpongeCPP`
- workflow filename: `publish-pypi.yml`
- environment: `pypi`

其中 `environment: pypi` 对 PyPI 来说不是强制项，但仓库当前 workflow
就是按这个名字设计的，建议保持一致。

## 发布步骤

1. 修改 `pyproject.toml` 中的版本号。
2. 先在本地做一次打包验证：

   ```bash
   pixi run python scripts/build_pypi.py
   ```

3. 提交并推送版本更新和相关 release note。
4. 在 GitHub 上创建对应版本的 Release。
5. 发布这个 GitHub Release。
6. GitHub Actions 会自动：
   - 构建所有配置好的 wheel
   - 构建 sdist
   - 汇总 artifacts
   - 通过 OIDC Trusted Publishing 发布到 PyPI

## 为什么暂时不切到 cibuildwheel

`cibuildwheel` 以后仍然值得考虑，但当前仓库先保留手写 workflow，
原因是它更容易审计和维护：

- `XpongeCPP` / `Xponge` 双包布局
- 最小 smoke 与完整验证的分层
- 当前显式的平台/架构矩阵

等发布矩阵和发布策略进一步稳定后，再考虑迁移到 `cibuildwheel` 会更合适。
