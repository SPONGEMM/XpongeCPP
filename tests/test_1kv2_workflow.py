from pathlib import Path

import XpongeCPP as Xponge


DATA_DIR = Path("/media/ylj/62dc0c74-e929-4dc8-8db9-632cb94b0cb8/Mokda_demos/1KV2/data")
PDB_1KV2_H = DATA_DIR / "1KV2_H.pdb"


def test_forcefield_imports_register_ff14sb_tip3p_and_ions():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    assert Xponge.has_template("MET")
    assert Xponge.has_template("SER")
    assert Xponge.has_template("WAT")
    assert Xponge.has_template("NA")
    assert Xponge.has_template("CL")
    assert Xponge.template_atom_count("WAT") == 3
    assert Xponge.template_atom_count("NA") == 1


def test_load_1kv2_h_pdb_preserves_protein_structure():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.load_pdb(str(PDB_1KV2_H))

    assert mol.atom_count == 5793
    assert mol.residue_count == 360
    assert mol.residues[0].name == "NMET"
    assert mol.residues[-1].name == "CSER"
    assert mol.residues[-1].name2atom("OXT").name == "OXT"
    assert mol.validate()


def test_1kv2_tip3p_ion_export_workflow(tmp_path):
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    import XpongeCPP.forcefield.amber.tip3p  # noqa: F401

    mol = Xponge.load_pdb(str(PDB_1KV2_H))
    water = Xponge.get_template_molecule("WAT")

    Xponge.Add_Solvent_Box(mol, water, 8.0, tolerance=2.5, n_solvent=64)
    Xponge.Add_Ions(mol, {"NA": 4, "CL": 2})

    counts = mol.residue_counts()
    assert counts["WAT"] == 12263
    assert counts["NA"] == 4
    assert counts["CL"] == 2
    assert mol.validate()

    outputs = Xponge.Save_SPONGE_Input(mol, prefix="spg", dirname=str(tmp_path))
    expected = {
        "LJ",
        "angle",
        "atom_name",
        "atom_type_name",
        "bond",
        "charge",
        "coordinate",
        "dihedral",
        "exclude",
        "mass",
        "nb14",
        "residue",
        "resname",
    }
    assert set(outputs) == expected

    atom_count = mol.atom_count
    residue_count = mol.residue_count
    assert (tmp_path / "spg_residue.txt").read_text().splitlines()[0] == f"{atom_count} {residue_count}"
    assert int((tmp_path / "spg_charge.txt").read_text().splitlines()[0]) == atom_count
    assert int((tmp_path / "spg_mass.txt").read_text().splitlines()[0]) == atom_count

    for key in ["bond", "angle", "dihedral", "nb14"]:
        lines = (tmp_path / f"spg_{key}.txt").read_text().splitlines()
        for line in lines[1:200]:
            fields = line.split()
            if not fields:
                continue
            width = {"bond": 2, "angle": 3, "dihedral": 4, "nb14": 2}[key]
            for item in fields[:width]:
                assert 0 <= int(item) < atom_count

    pdb_path = tmp_path / "spg.pdb"
    Xponge.save_pdb(mol, str(pdb_path))
    reloaded = Xponge.load_pdb(str(pdb_path))
    assert reloaded.atom_count == mol.atom_count
    assert reloaded.residue_count == mol.residue_count


def test_cpp_sources_are_split_by_subsystem():
    root = Path(__file__).resolve().parents[1] / "cpp"
    for subdir in ["core", "forcefield", "io", "topology", "solvation", "python"]:
        assert (root / subdir).is_dir(), subdir
        assert any((root / subdir).glob("*.cpp")), subdir
