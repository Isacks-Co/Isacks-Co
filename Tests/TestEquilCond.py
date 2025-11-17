# tests/test_equilibrium_condition.py
import sys
from Tests.TestBase import TestBase
import numpy as np

sys.path.append("..")
from SourceCode.Utils.equilibriumCondition import EquilibriumCondition


class TestequilibriumCondition(TestBase):

    # -------- checkInternalPressureStable --------

    def test_internal_pressure_pass(self):
        # Small oscillation + mean ~ 0  -> True
        list_pi = [-5e-4, 3e-4, 1e-4, -2e-4]
        self.assertTrue(EquilibriumCondition.checkInternalPressureStable(list_pi, tol=1e-3))

    def test_internal_pressure_fail_range(self):
        # Range too large -> False
        list_pi = [0.0, 0.0, 2.1e-3]
        self.assertFalse(EquilibriumCondition.checkInternalPressureStable(list_pi, tol=1e-3))

    def test_internal_pressure_fail_mean(self):
        # Mean far from 0 despite small range -> False
        list_pi = [0.0020, 0.0021, 0.0019]  # mean ≈ 2e-3 > tol_mean (1e-3)
        self.assertFalse(EquilibriumCondition.checkInternalPressureStable(list_pi, tol=1e-3))

    def test_internal_pressure_boundary_equal_tol(self):
        # Range equals tol AND mean within tol -> True
        list_pi = [0.0, 1e-3]  # delta == tol
        self.assertTrue(EquilibriumCondition.checkInternalPressureStable(list_pi, tol=1e-3))

    # -------- checkStable (generic) --------

    def test_stable_pass_absolute(self):
        # Small absolute oscillation -> True (provide a looser tol than default)
        xs = [1.0002, 0.9998, 1.0001]
        self.assertTrue(EquilibriumCondition.checkStable(xs, tol=1e-3))

    def test_stable_fail_absolute(self):
        xs = [0.0, 2e-5]
        self.assertFalse(EquilibriumCondition.checkStable(xs, tol=1e-5))

    def test_stable_boundary_equal_tol(self):
        xs = [0.0, 1e-5]
        self.assertTrue(EquilibriumCondition.checkStable(xs, tol=1e-5))

    def test_stable_single_value(self):
        # With your current implementation, a single value yields delta=0 -> True
        xs = [42.0]
        self.assertTrue(EquilibriumCondition.checkStable(xs, tol=1e-9))

    def test_stable_with_numpy_input(self):
        # Accepts numpy arrays too
        xs = np.array([10.0, 10.0005, 9.9996])
        self.assertTrue(EquilibriumCondition.checkStable(xs, tol=1e-3))


if __name__ == "__main__":
    import unittest
    unittest.main()
