"""Register Amber OL3 RNA templates."""

from ... import register_amber_parmdat_file, register_residue_templates_from_mol2_file
from . import data_path
from ._terminal_monophosphate import register_terminal_monophosphate

register_amber_parmdat_file(str(data_path("parm10.dat")))
register_residue_templates_from_mol2_file(str(data_path("RNA.mol2")))

from ._nucleic import configure_nucleic_templates

configure_nucleic_templates(("A", "C", "G", "U"))

register_terminal_monophosphate("rna", ("A", "C", "G", "U"))

CITATIONS = ("10.1021/ct200162x",)
