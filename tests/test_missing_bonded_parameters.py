from io import StringIO
import pytest
import XpongeCPP as x


def molecule(n,label):
    atoms='\n'.join(f'{i+1} C{i} {i*1.5} 0 0 {label} 1 ERR 0' for i in range(n))
    bonds='\n'.join(f'{i+1} {i+1} {i+2} 1' for i in range(n-1))
    return x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\nERR\n{n} {n-1} 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n{bonds}\n'))


def test_missing_bond_parameters_fail_instead_of_using_current_distance(tmp_path):
    m=molecule(2,'MISSINGBOND')
    with pytest.raises(RuntimeError,match='missing bond parameter'):
        x.save_sponge_input(m,prefix='input',dirname=str(tmp_path))


def test_missing_angle_parameters_fail_instead_of_using_generic_values(tmp_path):
    x.register_amber_bond_parameter('MISSINGANGLE','MISSINGANGLE',100.,1.5)
    m=molecule(3,'MISSINGANGLE')
    with pytest.raises(RuntimeError,match='missing angle parameter'):
        x.save_sponge_input(m,prefix='input',dirname=str(tmp_path))
