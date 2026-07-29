"""Reusable PACKMOL-Memgen lipid extension registration."""

from ... import (
    register_amber_frcmod_file,
    register_amber_parmdat_file,
    register_residue_templates_from_mol2_file,
)
from ...helper import Xprint
from . import data_path
from ._forcefield_family import require_forcefield_family
from ._lipid_common import configure_manifest


LIPID_EXTENSION_REFERENCE_TEXT = """Reference for Xponge lipid extensions:
1. PI, phosphoinositide, and LysoPL extensions
  Schott-Verdugo, S.; Gohlke, H.
    PACKMOL-Memgen: A Simple-To-Use, Generalized Workflow for Membrane-Protein-Lipid-Bilayer System Building.
    Journal of Chemical Information and Modeling 2019 59(6), 2522-2528.
    DOI: 10.1021/acs.jcim.9b00269

2. Supporting parameter families
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
"""


_REGISTERED = False


def register_lipid_extension():
    """Register the shared extension once after a supported lipid base."""
    global _REGISTERED
    require_forcefield_family("lipid", {"lipid17", "lipid21"})
    if _REGISTERED:
        return None
    register_amber_parmdat_file(str(data_path("glycam_06j/GLYCAM_06j.dat")))
    register_amber_frcmod_file(str(data_path("frcmod.lipid_ext")))
    register_residue_templates_from_mol2_file(str(data_path("lipid_ext.mol2")))
    manifest = configure_manifest(data_path("lipid_ext_manifest.json"))
    Xprint(LIPID_EXTENSION_REFERENCE_TEXT)
    _REGISTERED = True
    return manifest

