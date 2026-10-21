"""Soft LJ export must not request ordinary-LJ-dependent NB14 loading."""
from io import StringIO
import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep

def chain(n):
    label='SOFTEXPORT'
    x.register_amber_lj_parameter(label,label,.1,1.8)
    x.register_amber_bond_parameter(label,label,100.,1.5)
    x.register_amber_angle_parameter([label]*3,50.,1.9)
    x.register_amber_proper_dihedral_parameter([label]*4,2,2.,0.,reset=True)
    atoms='\n'.join(f'{i+1} C{i} {i*1.5} {i%2} .2 {label} 1 LIG 0.1' for i in range(n))
    bonds='\n'.join(f'{i+1} {i+1} {i+2} 1' for i in range(n-1))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nREG\n{n} {n-1} 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n{bonds}\n'))

def test_empty_soft_nb14_is_not_exported(tmp_path):
    m=fep.merge_force_field(chain(2),chain(2),.5)
    x.Save_SPONGE_Input(m,prefix='TMP',dirname=str(tmp_path))
    assert (tmp_path/'TMP_LJ_soft_core.txt').exists()
    assert not (tmp_path/'TMP_LJ.txt').exists()
    assert not (tmp_path/'TMP_nb14.txt').exists()

def test_unconverted_nonempty_soft_nb14_exports_exact_explicit_coefficients(tmp_path):
    m=chain(4)
    x.Save_SPONGE_Input(m,prefix='before',dirname=str(tmp_path))
    original=(tmp_path/'before_nb14.txt').read_text()
    _,_,kl,ke=map(float,original.splitlines()[1].split())
    m.enable_lj_soft_core()
    x.Save_SPONGE_Input(m,prefix='TMP',dirname=str(tmp_path))
    assert not (tmp_path/'TMP_nb14.txt').exists()
    rows=(tmp_path/'TMP_nb14_extra.txt').read_text().splitlines()
    assert int(rows[0])==1
    i,j,a,b,kee=map(float,rows[1].split())
    assert (i,j)==(0.,3.)
    assert [a,b,kee]==pytest.approx([kl*.1*3.6**12,kl*.2*3.6**6,ke],rel=2e-6)
    m.enable_lj_soft_core(False)
    x.Save_SPONGE_Input(m,prefix='after',dirname=str(tmp_path))
    assert (tmp_path/'after_nb14.txt').read_text()==original

def test_soft_explicit_pairs_are_preserved(tmp_path):
    m=fep.merge_force_field(chain(4),chain(4),.25)
    x.Save_SPONGE_Input(m,prefix='TMP',dirname=str(tmp_path))
    assert not (tmp_path/'TMP_nb14.txt').exists()
    rows=(tmp_path/'TMP_nb14_extra.txt').read_text().splitlines()
    assert int(rows[0])==1
    i,j,a,b,kee=map(float,rows[1].split())
    assert (i,j)==(0.,3.) and a>0 and b>0 and kee>0

def test_ordinary_nb14_contract_is_unchanged(tmp_path):
    x.Save_SPONGE_Input(chain(4),prefix='TMP',dirname=str(tmp_path))
    assert (tmp_path/'TMP_LJ.txt').exists()
    rows=(tmp_path/'TMP_nb14.txt').read_text().splitlines()
    assert int(rows[0])==1 and rows[1].split()[:2]==['0','3']
