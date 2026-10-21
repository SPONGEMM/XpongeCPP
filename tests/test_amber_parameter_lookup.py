import unittest
import XpongeCPP as xp


class AmberParameterLookupTests(unittest.TestCase):
  def test_native_lookup_uses_live_wildcards_and_preserves_type_case(self):
    snapshot = xp._core._snapshot_forcefield_registries()
    try:
        xp.register_amber_improper_dihedral_parameter(["X", "X", "Qcenter", "Qoxygen"], 2, 10.5, 180.)
        assert xp.amber_parameter_exists("IMPROPER", ["Qalpha", "qn", "Qcenter", "Qoxygen"])
        assert not xp.amber_parameter_exists("IMPROPER", ["Qalpha", "qn", "qcenter", "Qoxygen"])
        xp.register_amber_bond_parameter("QA", "qb", 100., 1.5)
        assert xp.amber_parameter_exists("BOND", ["qb", "QA"])
        assert not xp.amber_parameter_exists("BOND", ["QA", "QB"])
        xp.register_amber_angle_parameter(["QA", "QB", "qc"], 30., 110.)
        assert xp.amber_parameter_exists("ANGLE", ["qc", "QB", "QA"])
        xp.register_amber_proper_dihedral_parameter(["X", "QA", "QB", "X"], 3, 1., 0.)
        assert xp.amber_parameter_exists("DIHE", ["qa", "QA", "QB", "qd"])
        with self.assertRaises(ValueError):
            xp.amber_parameter_exists("BOND", ["QA"])
    finally:
        xp._core._restore_forcefield_registries(snapshot)
    assert not xp.amber_parameter_exists("IMPROPER", ["Qalpha", "qn", "Qcenter", "Qoxygen"])


if __name__ == "__main__":
    unittest.main()
