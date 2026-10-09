"""Register terminal phosphate variants with their owning Amber force field."""
from ... import configure_residue_template_tail, register_residue_templates_from_mol2_file
from . import data_path


def register_terminal_monophosphate(family, bases):
    """Keep the HOP3 monoanion separate from the default 5′-OH mapping."""
    register_residue_templates_from_mol2_file(
        str(data_path("terminal_monophosphate_" + family + ".mol2")))
    for base in bases:
        configure_residue_template_tail(base + "5MP", "O3'", 1.5, "C3'")
