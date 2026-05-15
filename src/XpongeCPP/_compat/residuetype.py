"""Legacy ResidueType algebra and handle helpers."""

from .._core import Molecule, get_template_molecule, has_template


def coerce_template_like(value):
    """Convert a legacy template-like object into a one-residue molecule."""
    if isinstance(value, Molecule):
        return value.deepcopy()
    if hasattr(value, "name") and has_template(value.name):
        return get_template_molecule(str(value.name)).deepcopy()
    raise TypeError("legacy residue algebra expects template-like values with a registered name")


def add_template_like(left, right):
    """Implement legacy `ResidueType + ResidueType` style sequence building."""
    return coerce_template_like(left) + coerce_template_like(right)


def repeat_template_like(value, count):
    """Implement legacy `ResidueType * int` style sequence building."""
    return coerce_template_like(value) * int(count)
