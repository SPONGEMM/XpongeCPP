"""Legacy helper.namespace compatibility bridge."""

from .._compat.imports import (
    remove_real_global_variable,
    set_global_alternative_names,
    set_real_global_variable,
    source,
)

__all__ = [
    "remove_real_global_variable",
    "set_global_alternative_names",
    "set_real_global_variable",
    "source",
]
