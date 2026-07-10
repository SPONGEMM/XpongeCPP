import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "lipid21"
LIPID17_DATA = Path(__file__).resolve().parent / "data" / "lipid17"
AMBER_DATA = ROOT / "src" / "XpongeCPP" / "data" / "amber"
TEMPLATE_NAMES = [
    "AR", "CHL", "DHA", "LAL", "MY", "OL", "PA", "PC", "PE", "PGR",
    "PGS", "PH-", "PS", "SA", "SPM", "ST",
]


def _run(code, *args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [sys.executable, "-c", code, *map(str, args)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )


def test_lipid21_packaged_data_and_special_connection_metadata():
    manifest = json.loads((AMBER_DATA / "lipid21_manifest.json").read_text())
    assert manifest["format_version"] == 2
    assert manifest["source_sha256"] == hashlib.sha256((DATA / "lipid21.lib").read_bytes()).hexdigest()
    assert manifest["source_license"] == "Public Domain (AmberTools dat/leap)"
    assert [entry["template"] for entry in manifest["templates"]] == TEMPLATE_NAMES
    entries = {entry["template"]: entry for entry in manifest["templates"]}
    assert entries["PA"]["head_link_conditions"][1]["parameter_degrees"] == -120.0
    assert entries["PA"]["tail_link_conditions"][1]["parameter_degrees"] == 120.0
    assert entries["SA"]["head_next_atom"] == "C13"
    assert entries["SPM"]["head_next_atom"] == "N11"
    assert entries["SPM"]["tail_next_atom"] == "C2"


def test_lipid21_import_registers_base_and_extension_once():
    result = _run(
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid21\n"
        "from XpongeCPP.forcefield.amber._lipid_ext import register_lipid_extension\n"
        "assert all(X.has_template(name) for name in ('PGS','SA','SPM','PI3'))\n"
        "register_lipid_extension()\n"
    )
    output = result.stdout + result.stderr
    assert output.count("Reference for Lipid21:") == 1
    assert output.count("Reference for Xponge lipid extensions:") == 1
    assert "Reference for Lipid17:" not in output


def test_all_lipid21_and_extension_templates_export_with_complete_parameters(tmp_path):
    _run(
        "import json, sys\n"
        "from pathlib import Path\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid21\n"
        "base = json.loads(Path(sys.argv[1]).read_text())\n"
        "ext = json.loads(Path(sys.argv[2]).read_text())\n"
        "for entry in base['templates'] + ext['templates']:\n"
        "    molecule = X.get_template_molecule(entry['template'])\n"
        "    assert molecule.atom_count == entry['atom_count']\n"
        "    assert abs(sum(atom.charge for atom in molecule.atoms) - entry['total_charge']) < 1e-7\n"
        "    output = X.Save_SPONGE_Input(molecule, prefix=entry['template'].replace('-', 'minus'), dirname=sys.argv[3])\n"
        "    assert {'bond','angle','dihedral','nb14'}.issubset(output)\n",
        AMBER_DATA / "lipid21_manifest.json",
        AMBER_DATA / "lipid_ext_manifest.json",
        tmp_path,
    )


def test_lipid21_representative_pdbs_and_extension_links(tmp_path):
    _run(
        "import sys\n"
        "from pathlib import Path\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid21\n"
        "data, ext, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])\n"
        "expected = {'PSM': (127,['PA','SPM','SA']), 'SSM': (133,['ST','SPM','SA']), 'PGS': (123,['PA','PGS','PA'])}\n"
        "for name, (count, residues) in expected.items():\n"
        "    molecule = X.load_pdb(str(data / (name + '.pdb')))\n"
        "    assert molecule.atom_count == count\n"
        "    assert [residue.name for residue in molecule.residues] == residues\n"
        "    assert len(molecule.residue_links) == 2\n"
        "    X.Save_SPONGE_Input(molecule, prefix=name, dirname=str(out))\n"
        "molecule = X.load_mmcif(str(data / 'PSM.cif'))\n"
        "assert molecule.atom_count == 127\n"
        "assert [residue.name for residue in molecule.residues] == ['PA','SPM','SA']\n"
        "assert len(molecule.residue_links) == 2\n"
        "X.Save_SPONGE_Input(molecule, prefix='PSM_cif', dirname=str(out))\n"
        "molecule = X.load_pdb(str(ext / 'AHPI3.pdb'))\n"
        "assert molecule.atom_count == 146 and len(molecule.residue_links) == 2\n"
        "X.Save_SPONGE_Input(molecule, prefix='AHPI3', dirname=str(out))\n",
        DATA,
        LIPID17_DATA,
        tmp_path,
    )


def test_lipid21_preserves_torsion_specific_nb14_scaling(tmp_path):
    _run(
        "import sys\n"
        "from pathlib import Path\n"
        "import XpongeCPP as X\n"
        "import XpongeCPP.forcefield.amber.lipid21\n"
        "X.Save_SPONGE_Input(X.get_template_molecule('PA'), prefix='pa', dirname=sys.argv[1])\n"
        "rows = (Path(sys.argv[1]) / 'pa_nb14.txt').read_text().splitlines()[1:]\n"
        "scales = {tuple(row.split()[2:]) for row in rows}\n"
        "assert ('0.166667', '0.833333') in scales\n"
        "assert ('0.500000', '0.833333') in scales\n",
        tmp_path,
    )
