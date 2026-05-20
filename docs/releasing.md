# XpongeCPP release publishing

This document describes the repository's current PyPI publishing path.

## Release strategy

The repository now uses two separate GitHub Actions workflows:

- `build-packages.yml`
  - push / pull request validation
  - builds wheels across the supported platform matrix
  - builds an sdist on Linux x64
  - runs metadata checks and smoke tests
- `publish-pypi.yml`
  - runs when a GitHub Release is published
  - rebuilds all wheels and the sdist
  - uploads all release artifacts to PyPI through Trusted Publishing

## Current platform matrix

Release wheels are built for:

- Linux x64
- Linux arm64
- macOS Intel
- macOS arm64
- Windows x64

An sdist is built on Linux x64.

Linux wheels are built through `cibuildwheel`, so the published Linux artifacts
carry PyPI-compatible manylinux platform tags instead of raw `linux_*` tags.

## Smoke-test split

Each wheel build runs a minimal smoke test:

1. create a fresh virtual environment
2. install `numpy`
3. install the built wheel with `--no-deps`
4. verify:
   - `import XpongeCPP`
   - `import Xponge`

This keeps the wheel matrix broad without requiring the entire optional chemistry
stack on every runner.

The fuller dependency install test remains part of the regular packaging
validation workflow on Linux x64.

## PyPI Trusted Publishing setup

The release workflow is designed for PyPI Trusted Publishing with GitHub
Actions.

Configure the PyPI project to trust this GitHub Actions publisher:

- owner: `SPONGEMM`
- repository: `XpongeCPP`
- workflow filename: `publish-pypi.yml`
- environment: `pypi`

The `environment: pypi` setting is optional from PyPI's perspective, but it is
strongly recommended and is what this repository expects.

## Release procedure

1. Update the version in `pyproject.toml`.
2. Run a local packaging check:

   ```bash
   pixi run python scripts/build_pypi.py
   ```

3. Commit and push the version bump and any release notes.
4. Create a GitHub Release for the version tag.
5. Publish the GitHub Release.
6. GitHub Actions will:
   - build all configured wheels
   - build the sdist
   - collect artifacts
   - publish them to PyPI through OIDC Trusted Publishing

## Why not cibuildwheel yet

The repository now uses `cibuildwheel` for Linux wheel publication while
keeping the native hand-written flow for macOS and Windows. This hybrid setup
keeps Linux wheels PyPI-compatible while still making the overall workflow easy
to audit:

- the dual `XpongeCPP` / `Xponge` wheel layout
- the minimal-smoke versus fuller-validation split
- the current explicit platform matrix

If the matrix grows further, moving more of the release flow to
`cibuildwheel` can still be revisited later.
