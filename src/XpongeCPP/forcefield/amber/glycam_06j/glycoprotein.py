"""Register GLYCAM-06j glycoprotein templates.

`HYP/NHYP/CHYP` are provided by the Amber protein force fields (`ff14sb`/`ff19sb`).
This module only defines the GLYCAM glycoprotein bridge residues and their
N/C-terminal forms, so users should load the protein FF together with GLYCAM.
"""

from .... import (
    configure_residue_template_head,
    configure_residue_template_tail,
    register_pdb_residue_name_mapping,
    register_residue_templates_from_mol2_file,
)
from .. import data_path
from . import *  # noqa: F401,F403

register_residue_templates_from_mol2_file(str(data_path("glycam_06j/glycoprotein.mol2")))

for resname in ("OLP", "OLT", "NLN", "OLS"):
    configure_residue_template_head(resname, "N", 1.3, "CA")
    configure_residue_template_head("C" + resname, "N", 1.3, "CA")
    configure_residue_template_tail(resname, "C", 1.3, "CA")
    configure_residue_template_tail("N" + resname, "C", 1.3, "CA")
    register_pdb_residue_name_mapping("head", resname, "N" + resname)
    register_pdb_residue_name_mapping("tail", resname, "C" + resname)
