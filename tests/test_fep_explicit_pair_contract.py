"""FEP-only regression: endpoint functions, explicit pairs, and fail-fast limits."""
from io import StringIO
import math

import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep


@pytest.mark.parametrize('namespace',['XpongeCPP','Xponge'])
def test_native_fep_helpers_resolve_in_both_public_namespaces(namespace):
    import importlib
    module=importlib.import_module(namespace+'.forcefield.special.fep')
    a=chain(5)
    assert module.intramolecule_nb_to_nb14(a,0) is a


def chain(count=5, charge=.2, label='FEPREG', extra_bond=None):
    x.register_amber_lj_parameter(label, label, .1, 1.8)
    x.register_amber_bond_parameter(label, label, 100., 1.5)
    x.register_amber_angle_parameter([label]*3, 50., 1.9)
    x.register_amber_proper_dihedral_parameter([label]*4, 2, 2., 0., reset=True)
    atoms='\n'.join(f'{i+1} C{i} {i*1.5} {i%2} .2 {label} 1 LIG {charge}' for i in range(count))
    edges = [(i, i + 1) for i in range(1, count)]
    if extra_bond is not None:
        edges.append(extra_bond)
    bonds='\n'.join(f'{i} {j} {k} 1' for i, (j, k) in enumerate(edges, 1))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nREG\n{count} {len(edges)} 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n{bonds}\n'))


def save(molecule, directory, name):
    x.Save_SPONGE_Input(molecule, prefix=name, dirname=str(directory))
    return {p.name.removeprefix(name+'_'):p.read_text() for p in directory.glob(name+'_*.txt')}


def rows(files, name):
    return [list(map(float,line.split())) for line in files[name].splitlines()[1:] if line.strip()]


def test_explicit_pair_endpoint_and_charge_squared_interpolation(tmp_path):
    a=chain(4,.2); b=chain(4,.6)
    a.add_nb14_extra(0,2,10.,20.,.4)
    b.add_nb14_extra(0,2,30.,40.,.8)
    before_a,before_b=save(a,tmp_path,'a'),save(b,tmp_path,'b')
    stored={}
    for index,lam in enumerate((0.,.5,1.,.25,.5,0.,1.)):
        out=fep.merge_force_field(a,b,lam,{'charge':lam*lam})
        files=save(out,tmp_path,f'm{index}')
        row=next(r for r in rows(files,'nb14_extra.txt') if r[:2]==[0.,2.])
        assert row[2:4]==pytest.approx([(1-lam)*10+lam*30,(1-lam)*20+lam*40])
        q=out.residues[0].atoms[0].charge
        expected=(1-lam)*.4*.2**2+lam*.8*.6**2
        assert row[4]*q*q==pytest.approx(expected,rel=2e-6,abs=1e-12)
        # Standard 1-4 also converts to extra: no duplicate ordinary NB14.
        assert 'nb14.txt' not in files
        assert save(a,tmp_path,'a')==before_a and save(b,tmp_path,'b')==before_b
        if lam in stored: assert files==stored[lam]
        stored[lam]=files


def test_intramolecular_conversion_excludes_once_and_is_idempotent(tmp_path):
    a=chain(5)
    before=save(a,tmp_path,'before')
    assert 'nb14_extra.txt' not in before
    assert fep.intramolecule_nb_to_nb14(a,a.residues[0]) is a
    first=save(a,tmp_path,'first')
    pair=rows(first,'nb14_extra.txt')
    assert len(pair)==1 and pair[0][:2]==[0.,4.]
    epsilon=.1; r=3.6
    assert pair[0][2:]==pytest.approx([epsilon*r**12,2*epsilon*r**6,1.],rel=2e-6)
    exclusions=first['exclude.txt'].splitlines()[1:]
    assert 4 in list(map(int,exclusions[0].split()[1:]))
    fep.intramolecule_nb_to_nb14(a,[0,0])
    assert save(a,tmp_path,'second')==first


def test_core_pair_transition_categories_are_read_from_endpoints():
    a = chain(6)
    b = a.deepcopy()
    topology = b._get_explicit_topology()
    assert any(row[:2] == (1, 4) for row in topology["nb14s"])
    topology["nb14s"] = [row for row in topology["nb14s"] if row[:2] != (1, 4)]
    topology["nb14s"].append((0, 4, .5, .8))
    topology["exclusions"][0].extend([4, 5])
    b._set_explicit_topology(topology)
    assert fep.core_pair_transitions(a, b, 0) == [
        (0, 4, "ordinary", "native14"), (0, 5, "ordinary", "excluded"),
        (1, 4, "native14", "excluded"),
    ]
    assert fep.core_pair_transitions(a, a, 0) == []


def test_dual_assembly_preserves_each_endpoint_core_exclusions():
    a, b = chain(6), chain(6, extra_bond=(2, 5))
    left, right, _ = fep.merge_dual_topology(a, 0, b, {i: i for i in range(6)})
    assert fep.core_pair_transitions(left, right, 0) == [
        (0, 4, "ordinary", "excluded"), (0, 5, "ordinary", "native14"),
        (1, 4, "native14", "excluded"), (1, 5, "ordinary", "excluded"),
    ]


@pytest.mark.parametrize('selector',[-1,1,99])
def test_bad_residue_selection_rejected_without_input_mutation(tmp_path,selector):
    a=chain();before=save(a,tmp_path,'before')
    with pytest.raises(IndexError): fep.intramolecule_nb_to_nb14(a,selector)
    with pytest.raises(IndexError): fep.get_free_molecule(a,selector)
    assert save(a,tmp_path,'after')==before


@pytest.mark.parametrize('lam',[float('nan'),float('inf'),-.01,1.01])
def test_invalid_lambda_rejected(lam):
    a=chain()
    with pytest.raises(ValueError,match='lambda'): fep.merge_force_field(a,a,lam)
    with pytest.raises(ValueError,match='lambda'): fep.merge_force_field(a,a,.5,{'charge':lam})


def test_unsupported_modes_do_not_silently_change_the_model(tmp_path):
    a=chain()
    with pytest.raises(NotImplementedError): fep.merge_force_field(a,a,.5,intra_fep=True)
    with pytest.raises(NotImplementedError): fep.get_free_molecule(a,0,intra_fep=True)
    with pytest.raises(ValueError,match='specific_lambda'): fep.merge_force_field(a,a,.5,{'typo':.2})
    # Harmonic improper is different from supported Amber periodic improper.
    a.add_improper_dihedral(0,1,2,3,10.,0.)
    assert 'improper_dihedral.txt' in save(a,tmp_path,'ordinary')
    with pytest.raises(ValueError,match='unsupported force group'): fep.merge_force_field(a,a,.5)


def test_zero_charge_product_is_explicitly_blocked_not_regularized():
    a=chain(4,.2);b=chain(4,-.2)
    with pytest.raises(ValueError,match='zero charge product'): fep.merge_force_field(a,b,.5)
    assert a.residues[0].atoms[0].charge==.2
    assert b.residues[0].atoms[0].charge==-.2


def test_charge_export_does_not_round_compensation_charge_to_zero(tmp_path):
    a=chain(4,1e-9)
    files=save(fep.merge_force_field(a,a,.5),tmp_path,'small')
    charges=list(map(float,files['charge.txt'].splitlines()[1:]))
    assert charges==pytest.approx([1e-9*18.2223]*4,rel=1e-12,abs=0.)


def test_ordinary_charge_export_unchanged(tmp_path):
    a=chain(4,.2)
    assert save(a,tmp_path,'normal')['charge.txt'].splitlines()[1:] == ['3.644460']*4


def test_standard_nb14_uses_endpoint_charge_products(tmp_path):
    a=chain(4,.2);b=chain(4,.6)
    x.register_amber_nb14_scale('FEPREG','FEPREG',.4,.7)
    files=save(fep.merge_force_field(a,b,.5,{'charge':.25}),tmp_path,'mixed')
    pair=next(row for row in rows(files,'nb14_extra.txt') if row[:2]==[0.,3.])
    assert pair[4]*.3**2 == pytest.approx(.7*(.2**2+.6**2)/2,rel=2e-6)
    assert pair[2:4]==pytest.approx([.4*.1*3.6**12,.4*2*.1*3.6**6],rel=2e-6)


def test_standard_and_explicit_same_pair_are_added_not_overwritten(tmp_path):
    a=chain(4,.2);b=chain(4,.6)
    x.register_amber_nb14_scale('FEPREG','FEPREG',.4,.7)
    a.add_nb14_extra(0,3,10.,20.,.4);b.add_nb14_extra(0,3,30.,40.,.8)
    files=save(fep.merge_force_field(a,b,1.),tmp_path,'combined')
    pairs=rows(files,'nb14_extra.txt')
    assert len(pairs)==1
    assert pairs[0][2:]==pytest.approx([30.+.4*.1*3.6**12,40.+.4*2*.1*3.6**6,1.5],rel=2e-6)


def test_bonded_specific_lambdas_are_independent(tmp_path):
    a=chain(4,.2,'FEPREGA');b=chain(4,.6,'FEPREGB')
    x.register_amber_bond_parameter('FEPREGB','FEPREGB',200.,1.5)
    x.register_amber_angle_parameter(['FEPREGB']*3,90.,1.9)
    x.register_amber_proper_dihedral_parameter(['FEPREGB']*4,2,6.,0.,reset=True)
    x.register_amber_nb14_scale('FEPREGA','FEPREGA',.4,.7)
    x.register_amber_nb14_scale('FEPREGB','FEPREGB',.4,.7)
    files=save(fep.merge_force_field(a,b,.5,{'charge':.25,'bond_base':.2,'angle':.7,'dihedral':.8}),tmp_path,'specific')
    assert [r[2] for r in rows(files,'bond.txt')]==pytest.approx([120.]*3)
    assert [r[3] for r in rows(files,'angle.txt')]==pytest.approx([78.]*2)
    assert [r[5] for r in rows(files,'dihedral.txt')]==pytest.approx([5.2])
    pair=rows(files,'nb14_extra.txt')[0]
    assert pair[4]*.3**2==pytest.approx(.7*(.2*.2**2+.8*.6**2),rel=2e-6)


def test_decoupling_keeps_physical_bonded_terms_and_source_input(tmp_path):
    a=chain(4,0.)
    before=save(a,tmp_path,'before')
    free=fep.get_free_molecule(a,0)
    for i,lam in enumerate((0.,.5,1.)):
        files=save(fep.merge_force_field(a,free,lam),tmp_path,f'decouple{i}')
        for name in ('bond.txt','angle.txt','dihedral.txt'):
            assert files[name]==before[name]
    assert save(a,tmp_path,'after')==before


def test_dual_topology_rejects_unsupported_force_loss():
    a=chain(4,0.);b=chain(4,0.)
    a.add_nb14_extra(0,2,1.,2.,0.)
    with pytest.raises(ValueError,match='explicit pairs'): fep.merge_dual_topology(a,0,b,{i:i for i in range(4)})
    c=chain(4,0.);c.add_improper_dihedral(0,1,2,3,10.,0.)
    with pytest.raises(ValueError,match='unsupported force group'): fep.merge_dual_topology(c,0,b,{i:i for i in range(4)})


def test_native_branch_lj_weight_leaves_charge_and_ordinary_lj_unchanged(tmp_path):
    a, b = chain(6, .2), chain(6, .2)
    x.register_amber_nb14_scale("FEPREG", "FEPREG", .4, .7)
    left, right, _ = fep.merge_dual_topology(a, 0, b, {0: 0, 1: 1})
    baseline = save(fep.merge_force_field(left, right, .5), tmp_path, "baseline")
    settings = {"nb14_junction_lj_a": 1., "nb14_internal_lj_a": .75}
    changed = save(fep.merge_force_field(left, right, .5, settings), tmp_path, "changed")
    before = {tuple(row[:2]): row[2:] for row in rows(baseline, "nb14_extra.txt")}
    after = {tuple(row[:2]): row[2:] for row in rows(changed, "nb14_extra.txt")}
    assert before.keys() == after.keys()
    assert after[(0., 3.)][:2] == pytest.approx([2 * value for value in before[(0., 3.)][:2]])
    assert after[(2., 5.)][:2] == pytest.approx([1.5 * value for value in before[(2., 5.)][:2]], rel=2e-6)
    assert [after[k][2] for k in after] == [before[k][2] for k in after]
    for name in ("charge.txt", "LJ_soft_core.txt", "bond.txt", "angle.txt", "dihedral.txt", "exclude.txt"):
        assert changed[name] == baseline[name]
