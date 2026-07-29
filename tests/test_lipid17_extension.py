import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "lipid17"
AMBER_DATA = ROOT / "src" / "XpongeCPP" / "data" / "amber"
TEMPLATE_NAMES = [
    "2-A", "2-B", "2-C", "CAM", "CLI", "ERG", "H2A", "H2B", "H2C",
    "P2A", "P2B", "P2C", "P2D", "P2E", "P2F", "P3-", "P3A", "P3B",
    "P3C", "P3D", "P3E", "P3F", "P3H", "PC1", "PC2", "PE1", "PE2",
    "PG1", "PG2", "PH3", "PH4", "PH5", "PI", "PI3", "PI4", "PI5",
    "SIT", "STI",
]


def _run_python(code, *args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [sys.executable, "-c", code, *map(str, args)],
        text=True,
        capture_output=True,
        cwd=ROOT,
        env=env,
        check=True,
    )


def _pdb_as_mmcif(path):
    tags = [
        "group_PDB", "id", "type_symbol", "label_atom_id", "auth_atom_id",
        "label_comp_id", "auth_comp_id", "label_asym_id", "auth_asym_id",
        "label_seq_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_alt_id",
        "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
        "pdbx_PDB_model_num",
    ]
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip()
        residue = line[17:20].strip()
        sequence = int(line[22:26])
        element = line[76:78].strip() or next(char for char in atom if char.isalpha())
        rows.append(
            f"ATOM {int(line[6:11])} {element.upper()} {atom} {atom} {residue} {residue} "
            f"A A {sequence} {sequence} ? . {float(line[30:38])} {float(line[38:46])} "
            f"{float(line[46:54])} 1.0 0.0 1"
        )
    return "data_lipid17\nloop_\n" + "\n".join(f"_atom_site.{tag}" for tag in tags) + "\n" + "\n".join(rows) + "\n#\n"


def test_packaged_lipid17_extension_artifacts_and_manifest():
    for filename in [
        "lipid_ext.mol2",
        "lipid_ext_manifest.json",
        "frcmod.lipid_ext",
        "lipid_ext.LICENSE.txt",
    ]:
        assert (AMBER_DATA / filename).is_file()

    manifest = json.loads((AMBER_DATA / "lipid_ext_manifest.json").read_text())
    assert manifest["source_sha256"] == "9c42d98a7eaa77edc2650958e8272ce9339d43502899cb8df4cffba914ceed59"
    assert manifest["template_count"] == 38
    assert [entry["template"] for entry in manifest["templates"]] == TEMPLATE_NAMES


def test_lipid17_import_prints_extension_references_once_and_not_resp():
    result = _run_python("import XpongeCPP.forcefield.amber.lipid17\n")
    output = result.stdout + result.stderr
    assert output.count("Reference for Lipid17:") == 1
    assert output.count("Reference for Xponge lipid extensions:") == 1
    assert "Reference for resp.py:" not in output


def test_lipid17_loads_glycam_parameters_without_glycam_templates():
    result = _run_python(
        "import json\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid17\n"
        "print(json.dumps({name: X.has_template(name) for name in ('PI', 'ROH', 'OME')}))\n"
    )
    state = json.loads(result.stdout.splitlines()[-1])
    assert state == {"PI": True, "ROH": False, "OME": False}


def test_all_extension_templates_export_with_complete_parameters(tmp_path):
    _run_python(
        "import json, sys\n"
        "from pathlib import Path\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid17\n"
        "manifest = json.loads(Path(sys.argv[1]).read_text())\n"
        "for entry in manifest['templates']:\n"
        "    molecule = X.get_template_molecule(entry['template'])\n"
        "    assert molecule.atom_count == entry['atom_count']\n"
        "    assert abs(sum(atom.charge for atom in molecule.atoms) - entry['total_charge']) < 1e-7\n"
        "    output = X.Save_SPONGE_Input(molecule, prefix=entry['template'].replace('-', 'minus'), dirname=sys.argv[2])\n"
        "    assert {'bond', 'angle', 'dihedral', 'nb14'}.issubset(output)\n",
        AMBER_DATA / "lipid_ext_manifest.json",
        tmp_path,
    )


def test_representative_pdbs_load_with_expected_links_and_export(tmp_path):
    _run_python(
        "import sys\n"
        "from pathlib import Path\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid17\n"
        "data, outdir = Path(sys.argv[1]), sys.argv[2]\n"
        "expected = {'AHPC': (140, ['AR','PC','DHA']), 'AHPI': (143, ['AR','PI','DHA']), "
        "'AHPI3': (146, ['AR','PI3','DHA']), 'AHPI45H': (151, ['AR','H2C','DHA']), "
        "'AHPI345H': (155, ['AR','P3H','DHA'])}\n"
        "for name, (atom_count, residue_names) in expected.items():\n"
        "    molecule = X.load_pdb(str(data / (name + '.pdb')))\n"
        "    assert molecule.atom_count == atom_count\n"
        "    assert [residue.name for residue in molecule.residues] == residue_names\n"
        "    assert len(molecule.residue_links) == 2\n"
        "    output = X.Save_SPONGE_Input(molecule, prefix=name, dirname=outdir)\n"
        "    assert {'bond', 'angle', 'dihedral', 'nb14', 'coordinate'}.issubset(output)\n",
        DATA,
        tmp_path,
    )


def test_representative_mmcif_loads_with_fragment_links(tmp_path):
    cif_path = tmp_path / "AHPI3.cif"
    cif_path.write_text(_pdb_as_mmcif(DATA / "AHPI3.pdb"))
    _run_python(
        "import sys\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid17\n"
        "molecule = X.load_mmcif(sys.argv[1])\n"
        "assert [residue.name for residue in molecule.residues] == ['AR', 'PI3', 'DHA']\n"
        "assert molecule.atom_count == 146\n"
        "assert len(molecule.residue_links) == 2\n"
        "X.Save_SPONGE_Input(molecule, prefix='AHPI3_cif', dirname=sys.argv[2])\n",
        cif_path,
        tmp_path,
    )


def test_lysophospholipid_templates_have_exactly_one_linkable_end(tmp_path):
    _run_python(
        "import sys\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid17\n"
        "residue_type = X.ResidueType.get_type\n"
        "for first, second in [('PA','PC1'),('PC2','PA'),('PA','PE1'),('PE2','PA'),('PA','PG1'),('PG2','PA')]:\n"
        "    molecule = residue_type(first) + residue_type(second)\n"
        "    assert len(molecule.residue_links) == 1\n"
        "    X.Save_SPONGE_Input(molecule, prefix=first + '_' + second, dirname=sys.argv[1])\n",
        tmp_path,
    )


def test_extension_preserves_distinct_amber_and_glycam_nb14_scales(tmp_path):
    _run_python(
        "import sys\n"
        "from pathlib import Path\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid17\n"
        "X.Save_SPONGE_Input(X.get_template_molecule('PI3'), prefix='pi3', dirname=sys.argv[1])\n"
        "rows = (Path(sys.argv[1]) / 'pi3_nb14.txt').read_text().splitlines()[1:]\n"
        "scales = {tuple(line.split()[2:]) for line in rows}\n"
        "assert ('0.500000', '0.833333') in scales\n"
        "assert ('1.000000', '1.000000') in scales\n",
        tmp_path,
    )
