"""Dual-topology ghosts retain auxiliary bonded terms, not physical LJ/q."""
import pytest
import XpongeCPP as x
from XpongeCPP.forcefield.special import fep
from test_fep_soft_bond_partition import chain, rows


@pytest.mark.parametrize('lam', [0., .25, .5, .75, 1.])
def test_dummy_branch_bond_angle_torsion_are_retained(tmp_path, lam):
    a, b = chain(4), chain(4)
    for atom in a.atoms:
        atom.charge = .1
    for atom in b.atoms:
        atom.charge = -.2
    left, right, _ = fep.merge_dual_topology(a, 0, b, {0: 0, 1: 1})
    mixed = fep.merge_force_field(left, right, lam)
    x.Save_SPONGE_Input(mixed, prefix='TMP', dirname=str(tmp_path))
    assert not (tmp_path / 'TMP_bond_soft.txt').exists()
    assert len(rows(tmp_path, 'TMP_bond.txt')) == 5
    assert len(rows(tmp_path, 'TMP_angle.txt')) == 4
    assert len(rows(tmp_path, 'TMP_dihedral.txt')) == 2
    assert all(r[2:] == [100., 1.5] for r in rows(tmp_path, 'TMP_bond.txt'))
    # Both branch torsions retain their complete source coefficients at all lambda.
    assert all(r[5] == pytest.approx(2.) for r in rows(tmp_path, 'TMP_dihedral.txt'))
    expected = [.1 * (1-lam) - .2 * lam] * 2 + [.1 * (1-lam)] * 2 + [-.2 * lam] * 2
    assert [atom.charge for atom in mixed.atoms] == pytest.approx(expected)
    assert (tmp_path / 'TMP_LJ_soft_core.txt').is_file()
    assert left.atoms[4].type == 'ZERO_LJ_ATOM'
    assert right.atoms[2].type == 'ZERO_LJ_ATOM'

    # A-only and B-only scaffolds must never acquire union-graph cross terms.
    # Check serialized atom identities: export may reorder disconnected groups.
    names = (tmp_path / 'TMP_atom_name.txt').read_text().splitlines()[1:]
    a_only = {i for i, name in enumerate(names) if name in {'C2', 'C3'}}
    b_only = {i for i, name in enumerate(names) if name in {'C2R2', 'C3R2'}}
    for filename, width in [('bond', 2), ('angle', 3), ('dihedral', 4), ('nb14_extra', 2)]:
        serialized = rows(tmp_path, f'TMP_{filename}.txt')
        for row in serialized:
            ids = set(map(int, row[:width]))
            assert not (ids & a_only and ids & b_only)
        # Same atom tuple with different harmonic minima/Fourier terms is valid;
        # an identical serialized function must not be emitted twice.
        functions = [tuple(row[:width] + row[width:]) for row in serialized]
        assert len(functions) == len(set(functions))
    excluded = set()
    for i, line in enumerate((tmp_path / 'TMP_exclude.txt').read_text().splitlines()[1:]):
        count, *others = map(int, line.split())
        assert count == len(others)
        excluded.update(tuple(sorted((i, j))) for j in others)
    assert {tuple(sorted((i, j))) for i in a_only for j in b_only} <= excluded


def test_common_parameter_change_does_not_require_soft_bond(tmp_path):
    a = chain(4, 'COMMON_A', 1.5, 100.)
    b = chain(4, 'COMMON_B', 1.6, 200.)
    left, right, _ = fep.merge_dual_topology(a, 0, b, {i: i for i in range(4)})
    mixed = fep.merge_force_field(left, right, .25)
    x.Save_SPONGE_Input(mixed, prefix='TMP', dirname=str(tmp_path))
    assert not (tmp_path / 'TMP_bond_soft.txt').exists()
    assert sorted(r[2:] for r in rows(tmp_path, 'TMP_bond.txt')) == [[50., 1.6]] * 3 + [[75., 1.5]] * 3
