"""Bootstrap helpers for installing the legacy compatibility runtime."""

from __future__ import annotations

from .aliases import install_top_level_aliases
from .assign import install_legacy_assign_patches
from .molecule import install_molecule_io_methods
from .runtime import install_legacy_runtime_patches
from .symbols import install_template_globals
from .workflows import ensure_mindsponge_todo_support


def install_legacy_bootstrap(
    namespace=None,
    *,
    install_aliases=False,
    install_templates=False,
    template_names=None,
    overwrite=False,
):
    """Install the first-wave legacy runtime in one explicit step."""
    install_legacy_runtime_patches(namespace=namespace)
    install_legacy_assign_patches()
    install_molecule_io_methods()
    ensure_mindsponge_todo_support()
    if install_aliases and namespace is not None:
        install_top_level_aliases(namespace)
    if install_templates and namespace is not None:
        install_template_globals(namespace=namespace, template_names=template_names, overwrite=overwrite)
    return namespace


Install_Legacy_Bootstrap = install_legacy_bootstrap

