"""The optimized backend must preserve endpoint ownership and FEP energies."""
from io import StringIO
import math

import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep


def _molecule(label, values):
    x.register_amber_lj_parameter(label, label, .1, 1.8)
    x.register_amber_bond_parameter(label, label, 100. if label == 'ISOA' else 200.,
                                    1.3 if label == 'ISOA' else 1.7)
    x.register_amber_angle_parameter([label]*3, 50., 1.9)
    for i, (k, phase, n) in enumerate(values):
        x.register_amber_proper_dihedral_parameter([label]*4, n, k, phase, reset=i == 0)
    return x.load_mol2(StringIO(f'''@<TRIPOS>MOLECULE
{label}
4 3 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C0 0.0 0.0 0.0 {label} 1 LIG 0.0
2 C1 1.0 1.0 0.2 {label} 1 LIG 0.0
3 C2 2.0 0.0 0.4 {label} 1 LIG 0.0
4 C3 3.0 1.0 0.6 {label} 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
'''))


def _save(molecule, directory, prefix):
    x.Save_SPONGE_Input(molecule, prefix=prefix, dirname=str(directory))
    return {p.name.removeprefix(prefix+'_'):p.read_text() for p in directory.glob(prefix+'_*.txt')}


def _torsion_energy(files, phi):
    return sum(float(row[5])*(1+math.cos(int(row[4])*phi-float(row[6])))
               for row in (line.split() for line in files['dihedral.txt'].splitlines()[1:]) if row)


@pytest.mark.parametrize('b_values', [
    [(4., 0., 1), (7., math.pi, 2)],
    [(7., math.pi, 2), (4., 0., 1)],
    [(4., .7, 3)],
])
def test_merge_is_input_preserving_order_independent_and_energy_correct(tmp_path, b_values):
    a = _molecule('ISOA', [(2., 0., 1), (3., math.pi, 2)])
    b = _molecule('ISOB', b_values)
    before_a, before_b = _save(a,tmp_path,'a'), _save(b,tmp_path,'b')
    stored = {}
    for iteration, lam in enumerate((0., .5, 1., .25, 1., .5, 0.)):
        result = fep.Merge_Force_Field(a,b,lam,{'charge':lam**2})
        files = _save(result,tmp_path,f'lambda{iteration}')
        assert _save(a,tmp_path,'a') == before_a
        assert _save(b,tmp_path,'b') == before_b
        values = [_torsion_energy(files,phi) for phi in (-2.1,-.4,0.,.8,2.9)]
        expected = [(1-lam)*_torsion_energy(before_a,phi)+lam*_torsion_energy(before_b,phi)
                    for phi in (-2.1,-.4,0.,.8,2.9)]
        assert values == pytest.approx(expected,abs=2e-5)
        # Each differing bond length retains its endpoint polynomial.
        bonds = [list(map(float,line.split())) for line in files['bond.txt'].splitlines()[1:] if line.strip()]
        for distance in (1.,1.3,1.7,2.):
            energy = sum(row[2]*(distance-row[3])**2 for row in bonds if {int(row[0]),int(row[1])} == {0,1})
            assert energy == pytest.approx((1-lam)*100*(distance-1.3)**2+lam*200*(distance-1.7)**2)
        if lam in stored:
            assert files == stored[lam]
        stored[lam] = files


def test_deepcopy_has_independent_atoms_and_preserves_normal_serialization(tmp_path):
    original = _molecule('ISOA', [(2.,0.,1),(3.,math.pi,2)])
    duplicate = original.deepcopy()
    before = _save(original,tmp_path,'original')
    assert _save(duplicate,tmp_path,'duplicate') == before
    duplicate.residues[0].atoms[0].charge = .75
    assert original.residues[0].atoms[0].charge == 0.
    assert _save(original,tmp_path,'original') == before
