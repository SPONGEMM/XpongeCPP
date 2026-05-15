"""Legacy namespace symbol injection helpers."""

from __future__ import annotations

import inspect
import sys

from .. import get_template_molecule, registered_template_names


def install_template_globals(namespace=None, template_names=None, overwrite=False):
    """Inject registered template molecules into a Python namespace."""
    if namespace is None:
        frame = inspect.currentframe()
        assert frame is not None and frame.f_back is not None
        namespace = frame.f_back.f_globals
    names = registered_template_names() if template_names is None else list(template_names)
    for name in names:
        if not overwrite and name in namespace:
            continue
        namespace[name] = get_template_molecule(name)
    return namespace


def sync_template_module_globals(template_names=None, overwrite=True):
    """Mirror registered template globals into loaded top-level compatibility modules."""
    names = registered_template_names() if template_names is None else list(template_names)
    for module_name in ("XpongeCPP", "Xponge"):
        module = sys.modules.get(module_name)
        if module is None:
            continue
        install_template_globals(namespace=module.__dict__, template_names=names, overwrite=overwrite)
    main_module = sys.modules.get("__main__")
    if main_module is not None:
        install_template_globals(namespace=main_module.__dict__, template_names=names, overwrite=overwrite)
    return names
