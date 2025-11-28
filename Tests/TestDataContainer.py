# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import logging
import numpy as np
import sys
import unittest
sys.path.append("../SourceCode")
from ASEWrappers import DataTrajectory, Frame
from TestBase import TestBase
from ase.lattice.cubic import FaceCenteredCubic

logger = logging.getLogger(__name__)


class TestDataContainer(TestBase):

    def setUp(self):
        super().setUp()
        atoms = FaceCenteredCubic(size=(5, 5, 5), symbol="Cu", pbc=True)
        self.traj = DataTrajectory(atoms)
        for i in range(0, 10):
            self.traj.append(Frame(0, {"Prop1": i}))

    def testLength(self):
        assert len(self.traj) == 10

    def testGetItem(self):
        assert isinstance(self.traj[0], Frame)
        for i, frame in enumerate(self.traj):
            assert frame["Prop1"] == i


class TestFrame(TestBase):

    def test(self):
        d = {"Prop1": 1, "Prop2": 2}
        frame = Frame(0, d)
        assert frame["Prop1"] == 1 and frame["Prop2"] == 2


if __name__ == "__main__":
    unittest.main()
