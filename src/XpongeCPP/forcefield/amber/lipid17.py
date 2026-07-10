"""Register Amber Lipid17 templates and parameters."""

from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("lipid", "lipid17")

import json

from ... import (
    configure_residue_template_head,
    configure_residue_template_tail,
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
)
from ...helper import Xprint
from . import data_path

register_amber_parmdat_file(str(data_path("lipid17.dat")))
register_residue_templates_from_mol2_file(str(data_path("lipid17.mol2")))
register_amber_parmdat_file(str(data_path("glycam_06j/GLYCAM_06j.dat")))
register_amber_frcmod_file(str(data_path("frcmod.lipid_ext")))
register_residue_templates_from_mol2_file(str(data_path("lipid_ext.mol2")))

for resname in "LAL PA MY OL ST AR DHA".split():
    configure_residue_template_head(resname, "C12", 1.5, "C13")
    configure_residue_template_tail(resname, "C12", 1.5, "C13")

for resname in "PC PE PS PGR PH-".split():
    configure_residue_template_head(resname, "C11", 1.5, "O11")
    configure_residue_template_tail(resname, "C21", 1.5, "O21")

with open(data_path("lipid_ext_manifest.json"), encoding="utf-8") as manifest_file:
    lipid_ext_manifest = json.load(manifest_file)

for entry in lipid_ext_manifest["templates"]:
    if entry["head_atom"] is not None:
        configure_residue_template_head(
            entry["template"], entry["head_atom"], 1.5, entry["head_next_atom"]
        )
    if entry["tail_atom"] is not None:
        configure_residue_template_tail(
            entry["template"], entry["tail_atom"], 1.5, entry["tail_next_atom"]
        )

Xprint("""Reference for Lipid17 and Xponge Lipid17 extensions:
1. Lipid14 / Lipid17 base parameters
  Dickson, C.J., Madej, B.D., Skjevik, A.A., Betz, R.M., Teigen, K., Gould, I.R., Walker, R.C.
    Lipid14: The Amber Lipid Force Field.
    Journal of Chemical Theory and Computation 2014 10(2), 865-879,
    DOI: 10.1021/ct4010307

2. PI, phosphoinositide, and LysoPL extensions
  Schott-Verdugo, S.; Gohlke, H.
    PACKMOL-Memgen: A Simple-To-Use, Generalized Workflow for Membrane-Protein-Lipid-Bilayer System Building.
    Journal of Chemical Information and Modeling 2019 59(6), 2522-2528.
    DOI: 10.1021/acs.jcim.9b00269

3. Supporting parameter families
  Kirschner, K.N. et al.
    GLYCAM06: A generalizable biomolecular force field. Carbohydrates.
    Journal of Computational Chemistry 2008 29(4), 622-655.
    DOI: 10.1002/jcc.20820

  Homeyer, N.; Horn, A.H.C.; Lanig, H.; Sticht, H.
    Revised AMBER Parameters for Bioorganic Phosphates.
    Journal of Chemical Theory and Computation 2012 8(11), 4405-4412.
    DOI: 10.1021/ct300613v

  Selected extension terms are derived from GAFF2 as documented in the source frcmod and Amber Reference Manual.
  GAFF lineage: Wang, J. et al., Journal of Computational Chemistry 2004 25, 1157-1174.
    DOI: 10.1002/jcc.20035
""")
