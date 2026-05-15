"""Legacy package-name shim that forwards Xponge imports to XpongeCPP."""

import sys

import XpongeCPP as _XpongeCPP
from XpongeCPP._compat.aliases import LEGACY_TOP_LEVEL_ALIAS_SPECS, install_top_level_aliases

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

for _name in dir(_XpongeCPP):
    if not _name.startswith("_") and _name not in _SHIM_SUBPACKAGES:
        globals().setdefault(_name, getattr(_XpongeCPP, _name))

install_top_level_aliases(globals())

_main = sys.modules.get("__main__")
if _main is not None:
    for _name, _value in list(globals().items()):
        if _name.startswith("_"):
            continue
        _main.__dict__.setdefault(_name, _value)

__all__ = list(getattr(_XpongeCPP, "__all__", [])) + list(LEGACY_TOP_LEVEL_ALIAS_SPECS)
