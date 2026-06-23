"""Amber force-field compatibility namespace."""

from .. import package_data_path
from ... import (
    configure_residue_template_head,
    configure_residue_template_tail,
    register_amber_frcmod_file,
    register_amber_nb14_scale,
    register_pdb_residue_name_mapping,
    set_lj_combining_rule,
)


set_lj_combining_rule("lorentz_berthelot")


def data_path(filename):
    return package_data_path("amber", filename)


def load_parameters_from_frcmod(filename, prefix=True):
    set_lj_combining_rule("lorentz_berthelot")
    register_amber_nb14_scale("X", "X", 0.5, 0.833333)
    return register_amber_frcmod_file(str(filename))


def configure_proline_like_terminal_mapping(resname, cterm_name, nterm_name=None):
    """Configure PDB terminal mappings for proline-like Amber residues."""

    configure_residue_template_head(resname, "N", 1.3, "CA")
    configure_residue_template_head(cterm_name, "N", 1.3, "CA")
    configure_residue_template_tail(resname, "C", 1.3, "CA")
    if nterm_name:
        configure_residue_template_tail(nterm_name, "C", 1.3, "CA")
    register_pdb_residue_name_mapping("head", resname, nterm_name or resname)
    register_pdb_residue_name_mapping("tail", resname, cterm_name)


Load_Parameters_From_Frcmod = load_parameters_from_frcmod
