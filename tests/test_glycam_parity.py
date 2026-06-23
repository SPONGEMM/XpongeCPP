from io import StringIO
from pathlib import Path

import XpongeCPP as Xponge


def _residue_name(residue):
    return getattr(residue, "name", None) or getattr(getattr(residue, "type", None), "name", None)


def test_ff14sb_hyp_c_terminal_mapping():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    pdb = StringIO(
        """\
ATOM      1  N   HYP A   1      -3.801  -7.214   0.158  1.00  0.00           N
ATOM      2  CD  HYP A   1      -2.343  -7.214   0.158  1.00  0.00           C
ATOM      3  CG  HYP A   1      -2.000  -5.850   0.680  1.00  0.00           C
ATOM      4  OD1 HYP A   1      -1.704  -5.934   2.097  1.00  0.00           O
ATOM      5  CB  HYP A   1      -3.218  -4.992   0.433  1.00  0.00           C
ATOM      6  CA  HYP A   1      -4.374  -5.881   0.158  1.00  0.00           C
ATOM      7  C   HYP A   1      -5.218  -5.692  -1.095  1.00  0.00           C
ATOM      8  O   HYP A   1      -5.328  -6.599  -1.916  1.00  0.00           O
ATOM      9  OXT HYP A   1      -6.363  -3.988  -4.755  1.00  0.00           O
TER
"""
    )
    mol = Xponge.load_pdb(pdb, ignore_hydrogen=True)

    assert _residue_name(mol.residues[0]) == "CHYP"


def test_ff19sb_hyp_n_terminal_mapping():
    import XpongeCPP.forcefield.amber.ff19sb  # noqa: F401

    pdb = StringIO(
        """\
ATOM      1  N   HYP A   1      -1.058  -0.265  -0.002  1.00  0.00           N
ATOM      2  CD  HYP A   1      -0.090  -1.353   0.078  1.00  0.00           C
ATOM      3  CG  HYP A   1       1.159  -0.666   0.543  1.00  0.00           C
ATOM      4  OD1 HYP A   1       1.294  -0.827   1.926  1.00  0.00           O
ATOM      5  CB  HYP A   1       0.991   0.791   0.183  1.00  0.00           C
ATOM      6  CA  HYP A   1      -0.442   1.044  -0.106  1.00  0.00           C
ATOM      7  C   HYP A   1      -0.866   1.702  -1.412  1.00  0.00           C
ATOM      8  O   HYP A   1      -1.621   1.120  -2.187  1.00  0.00           O
ATOM      9  N   ALA A   2       0.000   2.500  -1.800  1.00  0.00           N
ATOM     10  CA  ALA A   2       0.800   3.200  -2.700  1.00  0.00           C
ATOM     11  C   ALA A   2       1.900   2.400  -3.300  1.00  0.00           C
ATOM     12  O   ALA A   2       2.900   2.900  -3.800  1.00  0.00           O
ATOM     13  CB  ALA A   2       0.000   4.100  -3.700  1.00  0.00           C
TER
"""
    )
    mol = Xponge.load_pdb(pdb, ignore_hydrogen=True)

    assert Xponge.has_template("NHYP")
    assert _residue_name(mol.residues[0]) == "NHYP"


def test_glycoprotein_core_template_registration():
    import XpongeCPP.forcefield.amber.ff19sb  # noqa: F401
    import XpongeCPP.forcefield.amber.glycam_06j.glycoprotein  # noqa: F401

    expected = ["NLN", "NNLN", "CNLN", "OLS", "NOLS", "COLS", "OLT", "NOLT", "COLT", "OLP", "NOLP", "COLP"]
    for resname in expected:
        assert Xponge.has_template(resname)
    assert Xponge.has_template("HYP")
    assert Xponge.has_template("NHYP")
    assert Xponge.has_template("CHYP")
    bridge = Xponge.ResidueType.get_type("NOLP") + Xponge.ResidueType.get_type("CNLN")
    assert len(bridge.residue_links) == 1


def test_glycam_functional_group_templates_and_linkability():
    import XpongeCPP.forcefield.amber.glycam_06j.d_pyranose  # noqa: F401

    for resname in ["MEX", "SO3", "TBT", "6GB"]:
        assert Xponge.has_template(resname)

    mex_chain = Xponge.ResidueType.get_type("TBT") + Xponge.ResidueType.get_type("6GB") + Xponge.ResidueType.get_type("MEX")
    so3_chain = Xponge.ResidueType.get_type("TBT") + Xponge.ResidueType.get_type("6GB") + Xponge.ResidueType.get_type("SO3")
    assert len(mex_chain.residue_links) == 2
    assert len(so3_chain.residue_links) == 2


def test_glycam_coverage_audit_classifies_extension_layers():
    from XpongeCPP.forcefield.amber.glycam_06j.audit import audit_glycam_coverage

    report = audit_glycam_coverage(repo_root=Path(__file__).resolve().parents[1])
    for resname in ["MEX", "SO3", "TBT"]:
        assert resname in report["covered"]
    assert "CA2" in report["covered_elsewhere"]
    for resname in ["0AE", "0AF", "0GL", "0TV"]:
        assert resname in report["missing_modified_monosaccharides"]
