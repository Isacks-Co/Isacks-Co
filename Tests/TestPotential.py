import sys
import numpy as np
sys.path.append("..")
from Tests.TestBase import TestBase
from asap3.Internal.BuiltinPotentials import EMT
from asap3.Internal.BuiltinPotentials import LennardJones
from SourceCode.ASEWrappers import LennardJonesPotential, EMTPotential

import unittest

import logging

logger = logging.getLogger(__name__)

class TestPotential(TestBase):
    """Tests for the class Potential class. """

    def test_EMT(self):
        pot = EMTPotential()
        assert str(pot) == "EMT"
        
        assert isinstance(pot.getASEPotentialCalculator(),EMT)

    def test_LJ(self):
        
        pot = LennardJonesPotential(atomic_numbers= 1, sigmas=[0.5],epsilons = [0.5])

        assert str(pot) == "Lennard Jones"
        
        assert isinstance(pot.getASEPotentialCalculator(),LennardJones)
        
        
if __name__ == "__main__":
    unittest.main()