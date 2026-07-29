"""Register GLYCAM-06j L-pyranose templates and linkage semantics."""

from itertools import product

from .... import has_template, register_residue_templates_from_mol2_file
from .. import data_path
from . import *  # noqa: F401,F403

register_residue_templates_from_mol2_file(str(data_path("glycam_06j/l_pyranose.mol2")))

for code, terminal, conformation in product("ABCDEFGHJKLMNOPQRSTUVWXYZ".lower(), "012346789PQRSTUVWXYZ", "AB"):
    resname = terminal + code + conformation
    if not has_template(resname):
        continue
    if terminal == "1":
        configure_glycam_head(resname, 1)
    elif code in "cpbj":
        configure_glycam_tail(resname, "C2", "C3")
    elif code == "s":
        configure_glycam_tail(resname, "C2", "C3")
    else:
        configure_glycam_tail(resname, "C1", "C2")
    if terminal in "2ZYXTSP" or (terminal == "R" and code not in "VYW"):
        configure_glycam_head(resname, 2)
    elif terminal in "3WVQR":
        configure_glycam_head(resname, 3)
    elif terminal in "4U":
        configure_glycam_head(resname, 4)
    elif terminal != "1":
        configure_glycam_head(resname, int(terminal))
