"""Disjoint assembly must preserve and extend materialized physical forces."""
from io import StringIO
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep

def chain(n,label):
    x.register_amber_lj_parameter(label,label,.1,1.8)
    x.register_amber_bond_parameter(label,label,100.,1.5)
    atoms='\n'.join(f'{i+1} {label}{i} {i*1.5} 0 0 {label} 1 {label} 0' for i in range(n))
    bonds='\n'.join(f'{i+1} {i+1} {i+2} 1' for i in range(n-1))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nREG\n{n} {n-1} 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n{bonds}\n'))

def bond_names(m,path):
    x.Save_SPONGE_Input(m,prefix='TMP',dirname=str(path))
    names=[a.name for a in m.atoms]
    pairs={tuple(sorted((names[int(v[0])],names[int(v[1])]))) for line in (path/'TMP_bond.txt').read_text().splitlines()[1:] if (v:=line.split()) and abs(float(v[2]))>1e-12}
    assert int((path/'TMP_exclude.txt').read_text().split()[0])==m.atom_count
    return pairs

def test_materialized_left_appends_all_ordinary_right_forces(tmp_path):
    a=chain(2,'APPLEFT');b=chain(2,'APPRIGHT')
    mixed=fep.merge_force_field(a,a,.5)
    before=bond_names(mixed,tmp_path/'before')
    mixed.add_molecule(b)
    assert mixed.has_topology_override
    assert bond_names(mixed,tmp_path/'after')==before|{('APPRIGHT0','APPRIGHT1')}

def test_ordinary_left_preserves_materialized_right(tmp_path):
    a=chain(2,'APPSECOND');b=chain(2,'APPFIRST')
    mixed=fep.merge_force_field(a,a,.5)
    b.add_molecule(mixed)
    assert b.has_topology_override
    assert bond_names(b,tmp_path/'after')=={('APPFIRST0','APPFIRST1'),('APPSECOND0','APPSECOND1')}

def test_ordinary_append_still_rebuilds_without_materialized_cache(tmp_path):
    a=chain(2,'APPCONTROL1');a.add_molecule(chain(2,'APPCONTROL2'))
    assert not a.has_topology_override
    assert len(bond_names(a,tmp_path/'ordinary'))==2
