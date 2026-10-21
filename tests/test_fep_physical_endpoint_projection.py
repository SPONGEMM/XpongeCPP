"""FEP union IDs must not reselect an already chosen physical improper."""
from io import StringIO
import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep

def star():
    label='ENDPROJECT'
    x.register_amber_lj_parameter(label,label,.1,1.8)
    x.register_amber_bond_parameter(label,label,100.,1.5)
    x.register_amber_angle_parameter([label]*3,50.,1.9)
    x.register_amber_improper_dihedral_parameter([label]*4,2,1.1,3.141592653589793)
    atoms='\n'.join(f'{i+1} C{i} {i} {i%2} .2 {label} 1 LIG 0' for i in range(4))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nSTAR\n4 3 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n'))

def dihedrals(m,d):
    x.Save_SPONGE_Input(m,prefix='TMP',dirname=str(d))
    rows=[]
    for line in (d/'TMP_dihedral.txt').read_text().splitlines()[1:]:
        v=list(map(float,line.split()))
        if abs(v[5])>1e-12:rows.append((tuple(map(int,v[:4])),*v[4:]))
    return rows

@pytest.mark.parametrize('mapping',[{0:0,1:3,2:1,3:2},{0:2,1:0,2:3,3:1},{0:3,1:2,2:1,3:0}])
def test_B_endpoint_preserves_source_ordered_improper(tmp_path,mapping):
    a,b=star(),star();source=dihedrals(b,tmp_path/'physical')
    assert len(source)==1
    left,right,_=fep.merge_dual_topology(a,0,b,mapping)
    mixed=fep.merge_force_field(left,right,1.)
    actual=dihedrals(mixed,tmp_path/'dual')
    # Export groups hard-bond components; compare physical atom identities,
    # never raw pre-export union IDs against post-export row IDs.
    actual=[(tuple(mixed.atoms[i].name for i in ids),*values) for ids,*values in actual]
    expected=[(tuple(a.atoms[mapping[i]].name for i in ids),*values) for ids,*values in source]
    # Only whole-tuple reversal is equivalent, not arbitrary outer permutations.
    canonical=lambda rows:sorted((min(ids,ids[::-1]),*v) for ids,*v in rows)
    assert canonical(actual)==canonical(expected)

def test_ordinary_improper_is_unchanged_by_dual_construction(tmp_path):
    a,b=star(),star();before=dihedrals(b,tmp_path/'before')
    fep.merge_dual_topology(a,0,b,{0:0,1:3,2:1,3:2})
    assert dihedrals(b,tmp_path/'after')==before
