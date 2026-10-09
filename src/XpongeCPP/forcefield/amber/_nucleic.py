"""Shared Amber polymer boundaries; templates remain in the native registry."""

from ... import (
    configure_residue_template_head, configure_residue_template_tail,
    register_pdb_residue_name_mapping,
)


def configure_nucleic_templates(bases):
    """Configure internal, 5'-OH, 3'-OH and isolated nucleoside templates."""
    for base in bases:
        configure_residue_template_tail(base, "O3'", 1.5, "C3'")
        configure_residue_template_tail(base + "5", "O3'", 1.5, "C3'")
        configure_residue_template_head(base, "P", 1.5, "OP2")
        configure_residue_template_head(base + "3", "P", 1.5, "OP2")
        register_pdb_residue_name_mapping("head", base, base + "5")
        register_pdb_residue_name_mapping("tail", base, base + "3")
