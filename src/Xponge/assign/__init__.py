"""Legacy Xponge.assign shim."""

from XpongeCPP._compat.imports import extend_package_path

extend_package_path(globals(), "XpongeCPP.assign")

from XpongeCPP.assign import *  # noqa: F401,F403
