"""Shared residue-connection configuration for Amber lipid force fields."""

import json

from ... import configure_residue_template_head, configure_residue_template_tail


def configure_connection(residue_name, position, anchor, next_atom, length=1.5):
    """Configure one residue connection end."""
    configure = (
        configure_residue_template_head if position == "head"
        else configure_residue_template_tail
    )
    configure(residue_name, anchor, length, next_atom)


def configure_standard_chain(residue_name):
    for position in ("head", "tail"):
        configure_connection(residue_name, position, "C12", "C13")


def configure_standard_headgroup(residue_name):
    configure_connection(residue_name, "head", "C11", "O11")
    configure_connection(residue_name, "tail", "C21", "O21")


def configure_manifest(path):
    """Apply connection metadata from a generated lipid manifest."""
    with open(path, encoding="utf-8") as manifest_file:
        manifest = json.load(manifest_file)
    for entry in manifest["templates"]:
        for position in ("head", "tail"):
            anchor = entry[f"{position}_atom"]
            if anchor is not None:
                configure_connection(
                    entry["template"], position, anchor, entry[f"{position}_next_atom"]
                )
    return manifest

