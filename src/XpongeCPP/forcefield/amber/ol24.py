"""Amber OL24 DNA: parm10 plus the complete official OL24 correction."""

from ... import (
    register_amber_parmdat_file, register_amber_frcmod_file,
    register_residue_templates_from_mol2_file,
)
from . import data_path
from ._terminal_monophosphate import register_terminal_monophosphate
from ._nucleic import configure_nucleic_templates

register_amber_parmdat_file(str(data_path("parm10.dat")))
register_amber_frcmod_file(str(data_path("OL24.frcmod")))
register_residue_templates_from_mol2_file(str(data_path("ol24.mol2")))
configure_nucleic_templates(("DA", "DC", "DG", "DT"))

register_terminal_monophosphate("dna", ("DA", "DC", "DG", "DT"))

CITATIONS = ("10.1021/acs.jctc.4c01100",)
