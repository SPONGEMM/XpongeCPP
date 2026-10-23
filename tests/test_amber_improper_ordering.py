"""Named, parameter-slot AMBER improper ordering and term-wise FEP validation."""
from io import StringIO
import math
import numpy as np
import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep


@pytest.fixture(autouse=True)
def registry():
    saved=x._core._snapshot_forcefield_registries()
    yield
    x._core._restore_forcefield_registries(saved)


def setup(label):
    other=label+'3'
    for typ in (label,other):x.register_amber_lj_parameter(typ,typ,.1,1.8)
    for typ in (label,other):x.register_amber_bond_parameter(label,typ,100.,1.5)
    for typ in (label,other):x.register_amber_angle_parameter([label,label,typ],50.,1.9)
    return other


def star(label,order=(0,1,2,3),copies=1):
    other=label+'3';names=['C','A','B','D']
    points=[(0.,0.,0.),(1.,.2,.3),(0.,1.,0.),(-1.,-1.,.7)]
    atoms=[];bonds=[];ids={old:i+1 for i,old in enumerate(order)}
    for copy in range(copies):
        for old in order:
            xyz=points[old];typ=other if old==3 else label
            suffix=str(copy) if copies>1 else ''
            atoms.append(f'{4*copy+ids[old]} {names[old]+suffix} {xyz[0]+copy*10} {xyz[1]} {xyz[2]} {typ} 1 LIG 0')
        for old in range(1,4):bonds.append(f'{3*copy+old} {4*copy+ids[0]} {4*copy+ids[old]} 1')
    return x.load_mol2(StringIO('@<TRIPOS>MOLECULE\nSTAR\n'+f'{4*copies} {3*copies} 1\nSMALL\nUSER_CHARGES\n'
        '@<TRIPOS>ATOM\n'+'\n'.join(atoms)+'\n@<TRIPOS>BOND\n'+'\n'.join(bonds)+'\n'))


def terms(m,folder):
    x.Save_SPONGE_Input(m,prefix='TMP',dirname=str(folder))
    names=(folder/'TMP_atom_name.txt').read_text().splitlines()[1:]
    result=[]
    for line in (folder/'TMP_dihedral.txt').read_text().splitlines()[1:]:
        values=line.split();ids=tuple(names[int(i)] for i in values[:4])
        result.append((min(ids,ids[::-1]),int(values[4]),float(values[5]),float(values[6])))
    return sorted(result)


def test_parameter_slots_are_stable_under_atom_order_and_repeated_groups(tmp_path):
    label='DETCA';other=setup(label)
    x.register_amber_improper_dihedral_parameter([label,label,label,other],2,1.1,.37)
    expected=terms(star(label),tmp_path/'original')
    assert expected[0][0]==('A','B','C','D')
    assert terms(star(label,(3,2,0,1)),tmp_path/'permuted')==expected
    repeated=terms(star(label,copies=2),tmp_path/'copies')
    assert len(repeated)==2
    assert [(tuple(n[:-1] for n in ids),*v) for ids,*v in repeated]==[expected[0],expected[0]]


def test_multiple_periodicities_override_preserves_original_slots(tmp_path):
    label='DETMULTI';other=setup(label)
    x.register_amber_improper_dihedral_parameter([label,label,label,other],2,1.1,.37)
    x.register_amber_improper_dihedral_parameter([label,label,label,other],3,.4,.2)
    x.register_amber_improper_dihedral_parameter([other,label,label,label],2,2.2,.73)
    result=terms(star(label),tmp_path)
    assert sorted((n,k) for _,n,k,_ in result)==[(2,2.2),(3,.4)]
    assert next(ids for ids,n,*_ in result if n==2)==('B','C','A','D') # reverse(D,A,C,B)
    assert next(ids for ids,n,*_ in result if n==3)==('A','B','C','D')


def test_specificity_overrides_general_terms(tmp_path):
    label='DETGENERAL';other=setup(label)
    x.register_amber_improper_dihedral_parameter(['X','X',label,'X'],2,9.,.1)
    x.register_amber_improper_dihedral_parameter(['X','X',label,other],2,1.1,.2)
    result=terms(star(label,(3,2,0,1)),tmp_path)
    assert result==[(('A','B','C','D'),2,1.1,.2)]


def test_legacy_shared_registry_keeps_its_original_ordering(tmp_path):
    label='DETLEGACY';setup(label)
    x.register_amber_improper_dihedral_parameter(['X','X',label,'X'],2,1.1,.2,deterministic=False)
    assert terms(star(label),tmp_path)==[(('A','C','B','D'),2,1.1,.2)]


def test_opls_loader_explicitly_opts_out_of_amber_ordering(monkeypatch):
    from XpongeCPP.forcefield import opls
    output={key:'name\n' for key in ('LJ','bonds','angles','dihedrals','nb14')}
    output['periodic_impropers']='name phi0 k periodicity\nX-X-Q-X 180 4.184 2\n'
    calls=[]
    monkeypatch.setattr(opls,'load_ffitp',lambda _:output)
    monkeypatch.setattr(opls,'register_amber_improper_dihedral_parameter',lambda *a,**kw:calls.append((a,kw)))
    opls.load_parameter_from_ffitp('unused.itp','.')
    assert len(calls)==1 and calls[0][1]=={'deterministic':False}


def energy(records,xyz,lookup):
    value=0.
    for ids,n,k,phase in records:
        p=xyz[[lookup[name] for name in ids]]
        axis=p[2]-p[1];axis/=np.linalg.norm(axis)
        u=p[0]-p[1];v=p[3]-p[2]
        u-=np.dot(u,axis)*axis;v-=np.dot(v,axis)*axis
        phi=np.arctan2(np.dot(np.cross(axis,u),v),np.dot(u,v))
        value+=k*(1+np.cos(n*phi-phase))
    return value


@pytest.mark.parametrize('lam',[0.,.25,.5,.75,1.])
def test_fep_energy_force_preserves_physical_multiple_terms(tmp_path,lam):
    a_type='DETFEP_A';b_type='DETFEP_B';ao=setup(a_type);bo=setup(b_type)
    for n,k,p in [(2,1.1,.37),(3,.4,.2)]:
        x.register_amber_improper_dihedral_parameter([a_type,a_type,a_type,ao],n,k,p)
    # B intentionally has a different axis and non-special phases.
    for n,k,p in [(2,2.2,.73),(3,.8,.4)]:
        x.register_amber_improper_dihedral_parameter([bo,b_type,b_type,b_type],n,k,p)
    a=star(a_type);b=star(b_type,(3,2,0,1))
    original_a=terms(a,tmp_path/'a');original_b=terms(b,tmp_path/'b')
    amap={v.name:i for i,v in enumerate(a.atoms)}
    mapping={j:amap[v.name] for j,v in enumerate(b.atoms)}
    left,right,_=fep.merge_dual_topology(a,0,b,mapping)
    mixed=terms(fep.merge_force_field(left,right,lam),tmp_path/'mixed')
    xyz=np.array([[v.x,v.y,v.z] for v in a.atoms]);lookup={v.name:i for i,v in enumerate(a.atoms)}
    expected=lambda v:(1-lam)*energy(original_a,v,lookup)+lam*energy(original_b,v,lookup)
    assert energy(mixed,xyz,lookup)==pytest.approx(expected(xyz),abs=1e-6)
    for i in range(4):
        for j in range(3):
            l=xyz.copy();r=xyz.copy();l[i,j]-=1e-5;r[i,j]+=1e-5
            assert (energy(mixed,r,lookup)-energy(mixed,l,lookup))/2e-5==pytest.approx((expected(r)-expected(l))/2e-5,abs=1e-6)
    assert terms(a,tmp_path/'a_after')==original_a
    assert terms(b,tmp_path/'b_after')==original_b
