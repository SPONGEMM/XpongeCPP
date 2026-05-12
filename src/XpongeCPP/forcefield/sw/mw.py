"""Register packaged mW Stillinger-Weber templates and parameters."""

from ... import load_sw_parameter_file, register_template_molecule_from_mol2_file
from . import data_path

register_template_molecule_from_mol2_file(str(data_path("mw", "mW.mol2")))


def load_parameters(molecule, filename=None):
    if filename is None:
        filename = data_path("mw", "mW.sw")
    return load_sw_parameter_file(filename, molecule)
