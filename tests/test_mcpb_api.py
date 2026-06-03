from io import StringIO
import json
from pathlib import Path

import pytest
import XpongeCPP as Xponge
from XpongeCPP.mcpb import api as mcpb_api
from XpongeCPP.mcpb import seminario as mcpb_seminario
from XpongeCPP.qm.models import HessianResult


METAL_PDB_TEXT = """\
HETATM    1 ZN   MZN A   1       0.000   0.000   0.000  1.00  0.00          Zn
HETATM    2 N1   LIG A   2       2.100   0.000   0.000  1.00  0.00           N
END
"""

METAL_PDB_WITH_CONTEXT_TEXT = """\
HETATM    1 ZN   MZN A   1       0.000   0.000   0.000  1.00  0.00          Zn
HETATM    2 N1   LIG A   2       2.100   0.000   0.000  1.00  0.00           N
HETATM    3 O1   CTX A   3       5.000   0.000   0.000  1.00  0.00           O
END
"""

METAL_PDB_ANGLE_TEXT = """\
HETATM    1 ZN   MZN A   1       0.000   0.000   0.000  1.00  0.00          Zn
HETATM    2 N1   LIG A   2       2.100   0.000   0.000  1.00  0.00           N
HETATM    3 O1   OLG A   3       0.000   2.000   0.000  1.00  0.00           O
END
"""


def _load_metal_site():
    return Xponge.load_pdb(StringIO(METAL_PDB_TEXT))


def _load_metal_site_with_context():
    return Xponge.load_pdb(StringIO(METAL_PDB_WITH_CONTEXT_TEXT))


def _load_metal_site_with_angle():
    return Xponge.load_pdb(StringIO(METAL_PDB_ANGLE_TEXT))


def _load_single_donor_site(element: str):
    element = str(element)
    resname = f"M{element.upper():>2}"[-3:]
    text = (
        f"HETATM    1 {element.upper():>2}   {resname} A   1       0.000   0.000   0.000  1.00  0.00          {element:>2}\n"
        "HETATM    2 N1   LIG A   2       2.100   0.000   0.000  1.00  0.00           N\n"
        "END\n"
    )
    return Xponge.load_pdb(StringIO(text))


def test_mcpb_returns_structured_result_and_patches_residue_links():
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
    )

    assert isinstance(result, Xponge.MCPBResult)
    assert result.molecule is molecule
    assert result.request.ion_info[0].element == "Zn"
    assert result.selection.ion_atom_ids == (0,)
    assert result.selection.coordinating_atom_ids == (1,)
    assert result.selection.selected_residue_ids == (0, 1)
    assert result.connect_records == [(0, 1)]
    assert molecule.residue_links == [[0, 1]]
    assert "MZN" in result.registered_metal_templates
    assert Xponge.has_template("MZN")
    assert isinstance(result.small_model, Xponge.MCPBLocalModel)
    assert isinstance(result.large_model, Xponge.MCPBLocalModel)
    assert result.small_model.source_atom_ids == (0, 1)
    assert result.small_model.molecule.residue_count == 2
    assert result.small_model.molecule.residue_links == [[0, 1]]
    assert result.large_model.source_atom_ids == (0, 1)
    assert result.large_model.molecule.residue_count == 2
    assert result.sponge_ready is True
    assert any("provide basis" in item for item in result.pending_requirements)


def test_mcpb_requires_explicit_bonded_pairs_for_bonded_mode():
    molecule = _load_metal_site()

    with pytest.raises(ValueError, match="explicit bonded_pairs"):
        Xponge.MCPB(
            molecule,
            ion_ids=[0],
            ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        )


def test_mcpb_rejects_invalid_element_symbol():
    molecule = _load_metal_site()

    with pytest.raises(ValueError, match="unsupported element symbol"):
        Xponge.MCPB(
            molecule,
            ion_ids=[0],
            ion_info=[{"atom_id": 0, "element": "Xx", "formal_charge": 2, "spin": 1}],
            bonded_pairs=[(0, 1)],
        )


def test_mcpb_rejects_unsupported_explicit_ion_template_combination():
    molecule = _load_single_donor_site("Cu")

    with pytest.raises(ValueError, match="No packaged Amber ion template"):
        Xponge.MCPB(
            molecule,
            ion_ids=[0],
            ion_info=[{"atom_id": 0, "element": "Cu", "formal_charge": 3, "spin": 0}],
            bonded_pairs=[(0, 1)],
            method="blank",
        )


def test_mcpb_accepts_inferred_ion_info_but_warns_about_missing_charge_and_spin():
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        bonded_pairs=[(0, 1)],
    )

    assert result.request.ion_info[0].element == "Zn"
    assert any("formal_charge" in warning for warning in result.warnings)
    assert any("spin" in warning for warning in result.warnings)


def test_mcpb_large_model_can_include_additional_context_residues():
    molecule = _load_metal_site_with_context()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        additional_residue_ids=[2],
    )

    assert result.selection.selected_residue_ids == (0, 1, 2)
    assert result.small_model.source_atom_ids == (0, 1)
    assert result.small_model.molecule.residue_count == 2
    assert result.large_model.source_atom_ids == (0, 1, 2)
    assert result.large_model.molecule.residue_count == 3


def test_mcpb_can_patch_parent_charges_when_local_refit_runs(monkeypatch):
    molecule = _load_metal_site()

    def _fake_run_local_charge_refit(request, local_model):
        del request
        assert local_model.source_atom_ids == (0, 1)
        molecule.atoms[0].charge = 1.25
        molecule.atoms[1].charge = -1.25
        return {
            "charges": [1.25, -1.25],
            "updated_atom_ids": [0, 1],
            "assignment": None,
            "total_charge": 0,
            "spin": 1,
            "radius": {"Zn": 1.39},
        }

    monkeypatch.setattr(mcpb_api, "run_local_charge_refit", _fake_run_local_charge_refit)
    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        method="blank",
        basis="sto-3g",
    )

    assert result.updated_charge_atoms == [0, 1]
    assert molecule.atoms[0].charge == pytest.approx(1.25)
    assert molecule.atoms[1].charge == pytest.approx(-1.25)
    assert result.metadata["charge_refit"]["radius"] == {"Zn": 1.39}
    assert result.frcmod_path is not None


def test_mcpb_can_prepare_single_ion_export_for_sponge(tmp_path):
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
    )

    assert molecule.atoms[0].type == "Zn2+"
    assert molecule.atoms[0].mass == pytest.approx(65.4)
    assert result.sponge_ready is True
    Xponge.Save_SPONGE_Input(result.molecule, "case", dirname=str(tmp_path))
    assert (tmp_path / "case_LJ.txt").exists()
    assert "65.400" in (tmp_path / "case_mass.txt").read_text(encoding="utf-8")


def test_mcpb_blank_method_generates_frcmod_artifact():
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        method="blank",
    )

    assert result.frcmod_path is not None
    text = Path(result.frcmod_path).read_text(encoding="utf-8")
    assert "Xponge MCPB blank frcmod" in text
    assert "BOND" in text
    assert "N-Zn2+" in text or "Zn2+-N" in text


def test_mcpb_empirical_method_generates_nonzero_force_constants():
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        method="empirical",
    )

    assert result.frcmod_path is not None
    text = Path(result.frcmod_path).read_text(encoding="utf-8")
    assert "Xponge MCPB empirical frcmod" in text
    assert "BOND" in text
    assert "  0.0" not in text.split("ANGL")[0]


def test_mcpb_empirical_frcmod_is_applied_to_sponge_export(tmp_path):
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        method="empirical",
    )

    Xponge.Save_SPONGE_Input(result.molecule, "case", dirname=str(tmp_path))
    bond_text = (tmp_path / "case_bond.txt").read_text(encoding="utf-8")
    assert "46.600000" in bond_text


def test_mcpb_seminario_frcmod_is_applied_to_sponge_export(monkeypatch, tmp_path):
    molecule = _load_metal_site()

    def _fake_compute_hessian(assign, *, backend=None, basis="6-31g*", charge=0, spin=0, return_timings=False):
        del backend, basis, charge, spin, return_timings
        import numpy as np

        coords_bohr = np.asarray(assign.coordinates, dtype=float) / mcpb_seminario.BOHR_TO_ANGSTROM
        b_vector, _ = mcpb_seminario._bond_b_vector(coords_bohr, 0, 1)
        hessian_flat = 0.02 * np.outer(b_vector, b_vector)
        hessian = hessian_flat.reshape(2, 3, 2, 3).transpose(0, 2, 1, 3)
        return HessianResult(
            cartesian_hessian_au=hessian,
            coordinates_angstrom=[tuple(float(x) for x in row) for row in assign.coordinates],
            atom_symbols=list(assign.atoms),
            timings={"build": 0.0, "scf": 0.0, "hessian": 0.0, "total": 0.0},
        )

    monkeypatch.setattr(mcpb_seminario, "compute_hessian", _fake_compute_hessian)
    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        method="seminario",
        basis="sto-3g",
    )

    assert result.frcmod_path is not None
    frcmod_text = Path(result.frcmod_path).read_text(encoding="utf-8")
    assert "Xponge MCPB seminario frcmod" in frcmod_text
    Xponge.Save_SPONGE_Input(result.molecule, "case", dirname=str(tmp_path))
    bond_text = (tmp_path / "case_bond.txt").read_text(encoding="utf-8")
    assert "44.800000" in bond_text


def test_mcpb_seminario_can_register_angle_terms(monkeypatch, tmp_path):
    molecule = _load_metal_site_with_angle()

    def _fake_compute_hessian(assign, *, backend=None, basis="6-31g*", charge=0, spin=0, return_timings=False):
        del backend, basis, charge, spin, return_timings
        import numpy as np

        coords_bohr = np.asarray(assign.coordinates, dtype=float) / mcpb_seminario.BOHR_TO_ANGSTROM
        bond1, _ = mcpb_seminario._bond_b_vector(coords_bohr, 0, 1)
        bond2, _ = mcpb_seminario._bond_b_vector(coords_bohr, 0, 2)
        angle, _ = mcpb_seminario._angle_b_vector(coords_bohr, 1, 0, 2)
        hessian_flat = (
            0.02 * np.outer(bond1, bond1)
            + 0.03 * np.outer(bond2, bond2)
            + 0.10 * np.outer(angle, angle)
        )
        hessian = hessian_flat.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3)
        return HessianResult(
            cartesian_hessian_au=hessian,
            coordinates_angstrom=[tuple(float(x) for x in row) for row in assign.coordinates],
            atom_symbols=list(assign.atoms),
            timings={"build": 0.0, "scf": 0.0, "hessian": 0.0, "total": 0.0},
        )

    monkeypatch.setattr(mcpb_seminario, "compute_hessian", _fake_compute_hessian)
    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1), (0, 2)],
        method="seminario",
        basis="sto-3g",
    )

    assert result.frcmod_path is not None
    Xponge.Save_SPONGE_Input(result.molecule, "case", dirname=str(tmp_path))
    bond_text = (tmp_path / "case_bond.txt").read_text(encoding="utf-8")
    angle_text = (tmp_path / "case_angle.txt").read_text(encoding="utf-8")
    assert "48.700000" in bond_text
    assert "70.800000" in bond_text
    assert "91.840000" in angle_text


def test_mcpb_save_pdb_helper_writes_conect_records(tmp_path):
    molecule = _load_metal_site()
    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
    )

    outfile = tmp_path / "metal_site.pdb"
    Xponge.MCPB_Save_PDB(result.molecule, outfile)
    text = outfile.read_text(encoding="utf-8")
    assert "CONECT" in text


def test_mcpb_write_artifacts_helper_materializes_manifest_and_charge_overrides(tmp_path):
    molecule = _load_metal_site()

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 1}],
        bonded_pairs=[(0, 1)],
        method="blank",
        basis="sto-3g",
    )

    manifest = Xponge.MCPB_Write_Artifacts(result, tmp_path, prefix="metal_site")
    manifest_path = Path(manifest["manifest_path"])
    assert manifest_path.exists()
    saved_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert saved_manifest["pdb_path"].endswith("metal_site.pdb")
    assert saved_manifest["frcmod_path"].endswith("metal_site.frcmod")
    assert saved_manifest["charge_override_path"].endswith("metal_site.charge_overrides.json")
    assert "metal_site_LJ.txt" in saved_manifest["sponge_files"]
    charge_overrides = json.loads(Path(saved_manifest["charge_override_path"]).read_text(encoding="utf-8"))
    assert [record["atom_id"] for record in charge_overrides] == [0, 1]
    assert Path(saved_manifest["small_model_pdb"]).exists()
    assert Path(saved_manifest["large_model_pdb"]).exists()


@pytest.mark.parametrize(
    ("element", "formal_charge", "spin"),
    [
        ("Mg", 2, 0),
        ("Ca", 2, 0),
        ("Zn", 2, 0),
        ("Fe", 2, 0),
        ("Mn", 2, 0),
        ("Co", 2, 0),
        ("Ni", 2, 0),
        ("Cu", 2, 0),
    ],
)
def test_mcpb_tier1_metals_can_reach_blank_export_readiness(element, formal_charge, spin, tmp_path):
    molecule = _load_single_donor_site(element)

    result = Xponge.MCPB(
        molecule,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": element, "formal_charge": formal_charge, "spin": spin}],
        bonded_pairs=[(0, 1)],
        method="blank",
    )

    assert result.sponge_ready is True
    assert result.frcmod_path is not None
    assert molecule.atoms[0].type
    Xponge.Save_SPONGE_Input(result.molecule, "case", dirname=str(tmp_path))
    assert (tmp_path / "case_LJ.txt").exists()
