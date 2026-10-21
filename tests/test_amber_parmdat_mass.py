"""Amber atom masses must come from parameters, not element fallbacks."""
from io import StringIO
from pathlib import Path

import pytest
import XpongeCPP as xp


def fluorine_mol2():
    return StringIO("""@<TRIPOS>MOLECULE
MASS_F
1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 F1 0.0 0.0 0.0 f 1 MASS_F -0.25
@<TRIPOS>BOND
""")


def gaff_path():
    return Path(xp.__file__).parent/'data/amber/gaff.dat'


def test_parmdat_mass_reaches_mol2_and_template():
    xp.register_amber_parmdat_file(str(gaff_path()))
    molecule = xp.load_mol2(fluorine_mol2())
    assert molecule.atoms[0].mass == pytest.approx(19.0)
    assert xp.get_template_molecule('MASS_F').atoms[0].mass == pytest.approx(19.0)


def test_parmdat_mass_uses_supplied_value_and_last_registration(tmp_path):
    source = gaff_path().read_text()
    original = next(line for line in source.splitlines() if line.split()[:2] == ['f','19.00'])
    for mass in (17.125, 20.375):
        path = tmp_path/f'mass-{mass}.dat'
        path.write_text(source.replace(original, f'f {mass} 0.320 Fluorine', 1))
        xp.register_amber_parmdat_file(str(path))
        molecule = xp._core.load_mol2(fluorine_mol2())
        assert molecule.atoms[0].mass == pytest.approx(mass)


def test_frcmod_mass_override_remains_authoritative(tmp_path):
    xp.register_amber_parmdat_file(str(gaff_path()))
    path = tmp_path/'override.frcmod'
    path.write_text('mass override\nMASS\nf 21.750\n\n')
    xp.register_amber_frcmod_file(str(path))
    assert xp._core.load_mol2(fluorine_mol2()).atoms[0].mass == pytest.approx(21.750)
