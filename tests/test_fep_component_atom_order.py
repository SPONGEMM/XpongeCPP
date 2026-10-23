"""FEP serialization permutations must preserve every physical atom and term."""
from io import StringIO
import pytest
import XpongeCPP as x
from XpongeCPP import _core
from XpongeCPP.forcefield.special import fep

def star(label='FEPCOMPORDER',reverse=False):
    x.register_amber_lj_parameter(label,label,.1,1.8)
    x.register_amber_bond_parameter(label,label,100.,1.5)
    x.register_amber_angle_parameter([label]*3,50.,1.9)
    x.register_amber_improper_dihedral_parameter([label]*4,2,1.1,3.141592653589793)
    order=list(reversed(range(5))) if reverse else list(range(5))
    atoms='\n'.join(f'{i+1} C{i} {i*1.5} {i%2} .2 {label} 1 LIG 0.' for i in order)
    bonds='\n'.join(f'{i} 1 {i+1} 1' for i in range(1,5))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nREG\n5 4 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n{bonds}\n'))

@pytest.mark.parametrize('lam',[0.,.25,1.])
def test_fep_common_component_is_contiguous_with_explicit_identity_map(tmp_path,lam):
    a=star();b=star()
    left,right,_=fep.merge_dual_topology(a,0,b,{0:0,2:2,3:3,4:4})
    before=[(v.name,v.x,v.y,v.z) for v in left.atoms]
    mixed=fep.merge_force_field(left,right,lam)
    order=_core.prepare_sponge_atom_order(mixed)
    # Retained dummy bonds connect the complete scaffold: no artificial hard/
    # soft partition remains and the stable exporter must not reorder it.
    assert order==list(range(6))
    assert mixed.has_topology_override
    assert _core.prepare_sponge_atom_order(mixed)==list(range(6))
    x.Save_SPONGE_Input(mixed,prefix='TMP',dirname=str(tmp_path))
    assert [v.name for v in mixed.atoms]==['C0','C1','C2','C3','C4','C1R2']
    assert [(v.name,v.x,v.y,v.z) for v in left.atoms]==before
    rows=[list(map(float,r.split())) for r in (tmp_path/'TMP_bond.txt').read_text().splitlines()[1:]]
    assert rows==[[0.,1.,100.,1.5],[0.,2.,100.,1.5],[0.,3.,100.,1.5],
                  [0.,4.,100.,1.5],[0.,5.,100.,1.5]]
    assert not (tmp_path/'TMP_bond_soft.txt').exists()

def test_explicit_template_permutation_preserves_materialized_topology(tmp_path):
    a=star();b=star();mixed=fep.merge_force_field(a,b,.25)
    template=star(reverse=True)
    assert mixed.has_topology_override
    x.reorder_atoms_by_template(mixed,template)
    assert mixed.has_topology_override
    x.Save_SPONGE_Input(mixed,prefix='TMP',dirname=str(tmp_path))
    assert [v.name for v in mixed.atoms]==['C4','C3','C2','C1','C0']

def test_ordinary_non_fep_atom_order_unchanged(tmp_path):
    ordinary=star(reverse=True)
    before=[v.name for v in ordinary.atoms]
    assert _core.prepare_sponge_atom_order(ordinary)==list(range(5))
    x.Save_SPONGE_Input(ordinary,prefix='TMP',dirname=str(tmp_path))
    assert [v.name for v in ordinary.atoms]==before
