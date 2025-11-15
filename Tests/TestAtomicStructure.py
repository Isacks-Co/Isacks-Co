import sys
import numpy as np
sys.path.append("../SourceCode")
from Tests.TestBase import TestBase
from SourceCode.ASEWrappers import AtomicStructure
from SourceCode.ASEWrappers import EMTPotential

import unittest

import logging

logger = logging.getLogger(__name__)


class TestAtomicStructure(TestBase):
    """Tests for the class AtomicStructure. Plenty of room for more testing. Just some basic coverage right now."""
    def setUp(self):
        self.CU_struct = AtomicStructure.fromFile("TESTPOSCAR",pbc= True,supercells=[3,3,3],potential=EMTPotential())
        self.method_name = unittest.TestCase.id(self)
        logger.info(f"Running test: {self.method_name}")

    def testLength(self):
        assert len(self.CU_struct) == 27

    def testTemperature(self):
        assert self.CU_struct.temperature == 0
        self.CU_struct.setVelocitiesMB(300)
        assert np.isclose(self.CU_struct.temperature,300)
        self.CU_struct.setVelocitiesMB(0)
    
    def testVelocities(self):
        assert np.allclose(self.CU_struct.velocities,0)

    def testLabel(self):
        assert self.CU_struct.label[0:2] == "Cu"

    def testCohesive_energy(self):
        assert np.isclose(self.CU_struct.cohesive_energy,3.5,atol = 0.1 )
if __name__ == "__main__":
    unittest.main()
