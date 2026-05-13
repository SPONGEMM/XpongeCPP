import XpongeCPP as Xponge


def test_compat_module_adds_instance_style_molecule_save_methods(tmp_path):
    import XpongeCPP.compat  # noqa: F401
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401

    mol = Xponge.get_template_molecule("ALA")

    pdb_path = tmp_path / "lig.pdb"
    mol2_path = tmp_path / "lig.mol2"
    sponge_dir = tmp_path / "spg"
    sponge_dir.mkdir()

    mol.save_pdb(str(pdb_path))
    mol.save_mol2(str(mol2_path))
    outputs = mol.save_sponge_input(prefix="lig", dirname=str(sponge_dir))

    assert pdb_path.is_file()
    assert mol2_path.is_file()
    assert "bond" in outputs
    assert (sponge_dir / "lig_bond.txt").is_file()


def test_compat_layer_can_inject_legacy_template_globals():
    import XpongeCPP.forcefield.amber.ff14sb  # noqa: F401
    from XpongeCPP.compat import enable_legacy_namespace

    namespace = {}
    enable_legacy_namespace(namespace, template_names=["NALA", "ALA", "CALA"])

    mol = namespace["NALA"] + namespace["ALA"] * 12 + namespace["CALA"]

    assert mol.residue_count == 14
    assert [res.name for res in mol.residues[:2]] == ["NALA", "ALA"]
    assert mol.residues[-1].name == "CALA"
    assert mol.validate()
