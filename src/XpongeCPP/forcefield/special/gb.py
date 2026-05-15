"""Legacy generalized Born workflow shim."""

from ..._compat.workflows import Set_GB_Radius, set_gb_radius

SetGbRadius = set_gb_radius
setGbRadius = set_gb_radius

__all__ = ["set_gb_radius", "Set_GB_Radius", "SetGbRadius", "setGbRadius"]
