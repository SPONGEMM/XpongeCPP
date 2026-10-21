"""One-sided FEP bonds belong only to the soft table, not twice to energy."""
from io import StringIO
import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep

def chain(n,label='FEPSOFTPART',length=1.5,k=100.):
    x.register_amber_lj_parameter(label,label,.1,1.8)
    x.register_amber_bond_parameter(label,label,k,length)
    x.register_amber_angle_parameter([label]*3,50.,1.9)
    x.register_amber_proper_dihedral_parameter([label]*4,2,2.,0.,reset=True)
    atoms='\n'.join(f'{i+1} C{i} {i*1.5} {i%2} .2 {label} 1 LIG 0.' for i in range(n))
    bonds='\n'.join(f'{i+1} {i+1} {i+2} 1' for i in range(n-1))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nREG\n{n} {n-1} 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n{bonds}\n'))

def rows(folder,name):
    return [list(map(float,r.split())) for r in (folder/name).read_text().splitlines()[1:] if r.strip()]

@pytest.mark.parametrize('lam',[0.,.25,1.])
def test_one_sided_bonds_are_not_also_ordinary_bonds(tmp_path,lam):
    a=chain(2);b=chain(2)
    left,right,_=fep.merge_dual_topology(a,0,b,{0:0})
    mixed=fep.merge_force_field(left,right,lam)
    x.Save_SPONGE_Input(mixed,prefix='TMP',dirname=str(tmp_path))
    assert rows(tmp_path,'TMP_bond.txt')==[]
    soft=rows(tmp_path,'TMP_bond_soft.txt')
    assert sorted(soft)==[[0.,1.,100.,1.5,0.],[0.,2.,100.,1.5,1.]]

@pytest.mark.parametrize('lam',[0.,.25,1.])
def test_interleaved_dual_atoms_export_without_changing_atom_identity(tmp_path,lam):
    a=chain(4);b=chain(4)
    left,right,_=fep.merge_dual_topology(a,0,b,{0:0,2:2,3:3})
    names=[v.name for v in left.atoms]
    mixed=fep.merge_force_field(left,right,lam)
    x.Save_SPONGE_Input(mixed,prefix='TMP',dirname=str(tmp_path))
    assert [v.name for v in mixed.atoms]==names
    assert rows(tmp_path,'TMP_bond.txt')==[[2.,3.,100.,1.5]]
    soft=rows(tmp_path,'TMP_bond_soft.txt')
    assert len(soft)==4 and all(r[2:4]==[100.,1.5] for r in soft)
    assert [v.name for v in left.atoms]==names

def test_both_active_changed_length_stays_harmonic_and_lambda_weighted(tmp_path):
    a=chain(2,'FEPSOFTLEN_A',1.5,100.)
    b=chain(2,'FEPSOFTLEN_B',1.6,200.)
    mixed=fep.merge_force_field(a,b,.25)
    x.Save_SPONGE_Input(mixed,prefix='TMP',dirname=str(tmp_path))
    assert sorted(rows(tmp_path,'TMP_bond.txt'))==[[0.,1.,50.,1.6],[0.,1.,75.,1.5]]
    assert not (tmp_path/'TMP_bond_soft.txt').exists()
