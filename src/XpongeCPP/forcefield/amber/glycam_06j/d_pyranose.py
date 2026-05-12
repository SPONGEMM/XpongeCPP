"""Register GLYCAM-06j D-pyranose templates."""

from .... import register_residue_templates_from_mol2_file
from .. import data_path
from . import *  # noqa: F401,F403

register_residue_templates_from_mol2_file(str(data_path("glycam_06j/d_pyranose.mol2")))
