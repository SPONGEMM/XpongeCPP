"""Register GLYCAM-06j base carbohydrate templates and linkage semantics."""

from .... import (
    configure_residue_template_head,
    configure_residue_template_tail,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
)
from .. import data_path

register_amber_parmdat_file(str(data_path("glycam_06j/GLYCAM_06j.dat")))
register_residue_templates_from_mol2_file(str(data_path("glycam_06j/terminal.mol2")))
register_residue_templates_from_mol2_file(str(data_path("glycam_06j/functional_groups.mol2")))

configure_residue_template_head("ROH", "O1", 1.3, "HO1")
configure_residue_template_head("OME", "O", 1.3, "CH3")
configure_residue_template_head("MEX", "CH3", 1.52, "H1")
configure_residue_template_head("SO3", "S1", 1.74, "O1")
configure_residue_template_tail("TBT", "O1", 1.43, "C1")


def configure_glycam_head(resname, oxygen_index):
    """Set the glycosidic head atom to the indexed hydroxyl oxygen."""

    configure_residue_template_head(resname, f"O{oxygen_index}", 1.4, f"C{oxygen_index}")


def configure_glycam_tail(resname, atom_name, next_atom):
    """Set the glycosidic tail anchor for a carbohydrate residue."""

    configure_residue_template_tail(resname, atom_name, 1.4, next_atom)
