from pathlib import Path

import pytest

import XpongeCPP as Xponge


DATA_8RYK_DIR = Path(__file__).resolve().parent / "data" / "8ryk"


def _require_xpongelib():
    pytest.importorskip("XpongeLib")
    import XpongeCPP.forcefield.amber.gaff as gaff

    return gaff


def _residue_atom(residue, atom_name):
    name2atom = (
        getattr(residue, "name2atom", None)
        or getattr(residue, "Name2Atom", None)
        or getattr(residue, "Name2atom", None)
    )
    return name2atom(atom_name) if callable(name2atom) else getattr(residue, atom_name)


def _build_manual_8ryk_molecule():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    Xponge.load_mol2(DATA_8RYK_DIR / "edit_struct" / "CCS_3.gaff.mol2", as_template=True)
    mol = Xponge.Molecule("interactive_frcmod")

    phe_type = Xponge.ResidueType.get_type("PHE").deepcopy("PHE_1")
    phe_present = {
        "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2",
        "CZ", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HZ", "H",
    }
    phe_atoms = {
        str(getattr(atom, "name", "") or "").strip().upper()
        for atom in (getattr(Xponge.ResidueType.get_type("PHE"), "atoms", []) or [])
        if str(getattr(atom, "name", "") or "").strip()
    }
    phe_type.omit_atoms(sorted(phe_atoms - phe_present), charge=None)

    asp_type = Xponge.ResidueType.get_type("ASP").deepcopy("ASP_2")
    asp_present = {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "H", "HA", "HB2", "HB3"}
    asp_atoms = {
        str(getattr(atom, "name", "") or "").strip().upper()
        for atom in (getattr(Xponge.ResidueType.get_type("ASP"), "atoms", []) or [])
        if str(getattr(atom, "name", "") or "").strip()
    }
    asp_type.omit_atoms(sorted(asp_atoms - asp_present), charge=None)

    ccs_type = Xponge.ResidueType.get_type("CCS")

    gly_type = Xponge.ResidueType.get_type("GLY").deepcopy("GLY_4")
    gly_present = {"N", "CA", "C", "O", "H", "HA2", "HA3"}
    gly_atoms = {
        str(getattr(atom, "name", "") or "").strip().upper()
        for atom in (getattr(Xponge.ResidueType.get_type("GLY"), "atoms", []) or [])
        if str(getattr(atom, "name", "") or "").strip()
    }
    gly_type.omit_atoms(sorted(gly_atoms - gly_present), charge=None)

    mol.add_residue(Xponge.Residue(phe_type, directly_copy=True))
    mol.add_residue(Xponge.Residue(asp_type, directly_copy=True))
    mol.add_residue(Xponge.Residue(ccs_type, directly_copy=True))
    mol.add_residue(Xponge.Residue(gly_type, directly_copy=True))

    res1 = mol.residues[0]
    res2 = mol.residues[1]
    res3 = mol.residues[2]
    res4 = mol.residues[3]

    mol.add_residue_link(_residue_atom(res1, "N"), _residue_atom(res3, "CE"))
    mol.add_residue_link(_residue_atom(res2, "C"), _residue_atom(res3, "N"))
    mol.add_residue_link(_residue_atom(res3, "C"), _residue_atom(res4, "N"))
    return mol


def test_legacy_residuetype_deepcopy_and_omit_atoms_support_manual_script_pattern():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    phe_type = Xponge.ResidueType.get_type("PHE").deepcopy("PHE_NO_H")
    assert any(atom.name == "H" for atom in phe_type.atoms)

    phe_type.omit_atoms(["H"], charge=None)

    assert all(atom.name != "H" for atom in phe_type.atoms)

    mol = Xponge.Molecule("trimmed")
    mol.add_residue(Xponge.Residue(phe_type, directly_copy=True))

    assert mol.residue_count == 1
    assert mol.residues[0].name == "PHE_NO_H"
    assert all(atom.name != "H" for atom in mol.residues[0].atoms)


def test_parmchk2_gaff_accepts_molecule_objects(tmp_path):
    gaff = _require_xpongelib()
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    frcmod = tmp_path / "obj.frcmod"
    gaff.parmchk2_gaff(Xponge.get_template_molecule("PHE").deepcopy(), frcmod, direct_load=False, keep=True)

    assert frcmod.exists()
    text = frcmod.read_text()
    assert text.strip()
    assert "MASS" in text
    assert "NONBON" in text


def test_xpongecpp_matches_manual_8ryk_script_contract(tmp_path):
    gaff = _require_xpongelib()

    raw_mol2 = tmp_path / "manual.raw.mol2"
    raw_frcmod = tmp_path / "manual.raw.frcmod"
    mol = _build_manual_8ryk_molecule()

    Xponge.Save_Mol2(mol, raw_mol2)
    gaff.parmchk2_gaff(raw_mol2, raw_frcmod, direct_load=False, keep=True)

    assert raw_mol2.exists()
    assert raw_frcmod.exists()
    assert raw_frcmod.stat().st_size > 0
    frcmod_text = raw_frcmod.read_text()
    assert "MASS" in frcmod_text
    assert "NONBON" in frcmod_text
