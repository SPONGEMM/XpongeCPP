"""Register packaged silicon EDIP templates and parameters."""

from ... import load_edip_parameter_file, register_template_molecule_from_mol2_file
from . import data_path

register_template_molecule_from_mol2_file(str(data_path("si", "Si.mol2")))


def load_parameters(molecule, filename=None):
    if filename is None:
        filename = data_path("si", "Si.edip")
    return load_edip_parameter_file(filename, molecule)
