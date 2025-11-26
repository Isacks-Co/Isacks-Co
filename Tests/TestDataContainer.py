import sys
import numpy as np
sys.path.append("../SourceCode")
from TestBase import TestBase

from ase.lattice.cubic import FaceCenteredCubic

import unittest

import logging

logger = logging.getLogger(__name__)

from ASEWrappers import DataTrajectory,Frame

class TestDataContainer(TestBase):
    
    def setUp(self):
        super().setUp()
        atoms = FaceCenteredCubic(size=(5, 5, 5), symbol="Cu", pbc=True)
        self.traj = DataTrajectory(atoms)
        for i in range(0,10):
            self.traj.append(Frame(0, {"Prop1":i}))

    def testLength(self):
        assert len(self.traj) == 10

    def testGetItem(self):
        assert isinstance(self.traj[0],Frame) 
        for i,frame in enumerate(self.traj):
            assert frame["Prop1"] == i
    
    




class TestFrame(TestBase):
    
    def test(self):
        d = {"Prop1":1,"Prop2":2}
        frame = Frame(0, d)
        assert frame["Prop1"] == 1 and frame["Prop2"] == 2
    
    
if __name__ == "__main__":
    unittest.main()