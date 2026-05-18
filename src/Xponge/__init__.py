"""Legacy package-name shim that forwards Xponge imports to XpongeCPP."""

import XpongeCPP as _XpongeCPP
from XpongeCPP._compat.aliases import LEGACY_TOP_LEVEL_ALIAS_SPECS, install_top_level_aliases
from XpongeCPP._compat.imports import copy_public_attributes, install_main_namespace_exports

from XpongeCPP import *  # noqa: F401,F403

_SHIM_SUBPACKAGES = {
    "analysis",
    "assign",
    "build",
    "forcefield",
    "helper",
    "load",
    "mdrun",
    "process",
    "tools",
}

copy_public_attributes(_XpongeCPP, globals(), skip=_SHIM_SUBPACKAGES)

install_top_level_aliases(globals())
install_main_namespace_exports(globals())

__all__ = list(getattr(_XpongeCPP, "__all__", [])) + list(LEGACY_TOP_LEVEL_ALIAS_SPECS)
