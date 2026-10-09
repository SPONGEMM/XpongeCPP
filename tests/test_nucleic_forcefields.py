import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run(code, *args):
    result = subprocess.run([sys.executable, "-c", code, *map(str, args)],
                            cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("module,first,last,charge", [
    ("ol3", "A5", "U3", -1), ("bsc1", "DA5", "DT3", -1),
    ("ol15", "DA5", "DT3", -1), ("ol24", "DA5", "DT3", -1),
    ("ol3", "A5MP", "U3", -2), ("ol24", "DA5MP", "DT3", -2),
])
def test_complete_nucleic_chain_exports_with_integer_charge(tmp_path, module, first, last, charge):
    run("""
import importlib, sys
import XpongeCPP as X
importlib.import_module('XpongeCPP.forcefield.amber.' + sys.argv[1])
rt = X.ResidueType.get_type
m = rt(sys.argv[2]) + rt(sys.argv[3])
assert len(m.residue_links) == 1
assert abs(sum(a.charge for a in m.atoms) - float(sys.argv[4])) < 1e-6
X.Save_SPONGE_Input(m, prefix='nucleic', dirname=sys.argv[5], format='bundle')
import h5py
with h5py.File(sys.argv[5] + '/nucleic_topology.spgt.h5') as f:
    assert f['/topology/atom_count'][()] == len(m.atoms)
""", module, first, last, charge, tmp_path)


@pytest.mark.parametrize("modules", [("ol3", "ol24"), ("ol24", "ol3")])
def test_terminal_templates_follow_owning_forcefield_in_either_import_order(modules):
    run("""
import importlib, sys
import XpongeCPP as X
importlib.import_module('XpongeCPP.forcefield.amber.' + sys.argv[1])
get = X.get_template_molecule
first, absent = ('A5MP', 'DA5MP') if sys.argv[1] == 'ol3' else ('DA5MP', 'A5MP')
assert {'P', 'OP3', 'HOP3'} <= {a.name for a in get(first).atoms}
try:
    get(absent)
except (KeyError, ValueError, RuntimeError, IndexError):
    pass
else:
    raise AssertionError('loaded templates from an unrelated forcefield')
importlib.import_module('XpongeCPP.forcefield.amber.' + sys.argv[2])
for base in ('A', 'DA'):
    assert 'P' not in {a.name for a in get(base + '5').atoms}
    assert {'P', 'OP3', 'HOP3'} <= {a.name for a in get(base + '5MP').atoms}
""", *modules)

@pytest.mark.parametrize("water", ["tip3p", "spce", "tip4pew", "opc"])
@pytest.mark.parametrize("fmt", ["pdb", "mmcif"])
def test_single_atom_zinc_alias_keeps_source_atom_identity(tmp_path, water, fmt):
    run("""
import importlib,sys
from io import StringIO
import XpongeCPP as X
importlib.import_module('XpongeCPP.forcefield.amber.' + sys.argv[1])
if sys.argv[2] == 'pdb':
    text = 'HETATM   17 ZN    ZN A 101       0.000   0.000   0.000  1.00  0.00          ZN\\nEND\\n'
    m = X.load_pdb(StringIO(text), infer_terminals=False)
else:
    text = '''data_zinc
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 17 Zn ZN ZN A 101 0 0 0 101 ZN A ZN 1
#
'''
    m = X.load_mmcif(StringIO(text), infer_terminals=False)
a = m.atoms[0]
assert (a.name,a.serial,a.type,a.charge) == ('ZN',17,'Zn2+',2.0), (a.name,a.serial,a.type,a.charge)
X.Save_SPONGE_Input(m,prefix='zinc',dirname=sys.argv[3],format='bundle')
""", water, fmt, tmp_path)


def test_off_assets_and_charges_are_reproducible():
    spec = importlib.util.spec_from_file_location("converter", ROOT / "scripts/convert_amber_off_to_mol2.py")
    converter = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = converter
    spec.loader.exec_module(converter)
    data = ROOT / "src/XpongeCPP/data/amber"
    for source, target, suffix in [("ff-nucleic-OL24.lib", "ol24", ""), ("terminal_monophosphate.lib", "terminal_monophosphate", "MP")]:
        path = data / "amber_sources" / source
        templates = converter.parse_off(path)
        for template in templates:
            template.name += suffix
        assert converter.render_mol2(path, templates) == (data / (target + ".mol2")).read_text()
