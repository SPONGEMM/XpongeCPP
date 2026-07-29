"""Register Amber Lipid21 and the shared Xponge lipid extensions."""

from ... import register_amber_parmdat_file, register_residue_templates_from_mol2_file
from ...helper import Xprint
from . import data_path
from ._forcefield_family import activate_forcefield_family
from ._lipid_common import configure_manifest
from ._lipid_ext import register_lipid_extension


activate_forcefield_family("lipid", "lipid21")

register_amber_parmdat_file(str(data_path("lipid21.dat")))
register_residue_templates_from_mol2_file(str(data_path("lipid21.mol2")))
lipid21_manifest = configure_manifest(data_path("lipid21_manifest.json"))

Xprint("""Reference for Lipid21:
  Dickson, C.J.; Walker, R.C.; Gould, I.R.
    Lipid21: Complex Lipid Membrane Simulations with AMBER.
    Journal of Chemical Theory and Computation 2022 18(3), 1726-1736.
    DOI: 10.1021/acs.jctc.1c01217
""")

register_lipid_extension()
