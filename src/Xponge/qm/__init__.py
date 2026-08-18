"""Legacy Xponge.qm namespace backed by XpongeCPP."""

from XpongeCPP._compat.imports import extend_package_path

extend_package_path(globals(), "XpongeCPP.qm")

from XpongeCPP.qm import *  # noqa: F401,F403
