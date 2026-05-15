"""Legacy Xponge.helper shim."""

from XpongeCPP._compat.imports import extend_package_path

extend_package_path(globals(), "XpongeCPP.helper")

from XpongeCPP.helper import *  # noqa: F401,F403
