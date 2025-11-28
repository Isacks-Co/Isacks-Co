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
from ASEWrappers import Integrator, VelocityVerletIntegrator, LangevinIntegrator, IsotropicMTKNPTIntegrator, \
    BerendsenNPTIntegrator
from TestBase import TestBase


logger = logging.getLogger(__name__)


class TestIntegrator(TestBase):
    """Tests for the class Integrator class. """

    def generalTest(self):
        integrator = Integrator(1)
        assert len(integrator.attachments) == 0
        integrator.attach(1, 1)
        assert len(integrator.attachments) == 1

    def testVelocityVerlet(self):
        integrator = VelocityVerletIntegrator(1)
        assert integrator.ensemble == "NVE"

        assert str(integrator.integrator_partial).startswith("functools.partial(<function VelocityVerlet")

    def testLangevin(self):
        integrator = LangevinIntegrator(1, 300, 0.01)

        assert integrator.ensemble == "NVT"
        print(integrator.integrator_partial)
        assert str(integrator.integrator_partial).startswith("functools.partial(<function Langevin")

    def testIsotropicMTK(self):
        integrator = IsotropicMTKNPTIntegrator(1, 300, 0, 1, 1)
        assert integrator.ensemble == "NPT"
        print(integrator.integrator_partial)
        assert str(integrator.integrator_partial).startswith(
            "functools.partial(<class 'asap3.md.nose_hoover_chain.IsotropicMTKNPT'>")

    def testNPTBerendsen(self):
        integrator = BerendsenNPTIntegrator(1, 300, 0, 1e-7)
        assert integrator.ensemble == "NPT"
        print(integrator.integrator_partial)
        assert str(integrator.integrator_partial).startswith("functools.partial(<class 'asap3.md.nptberendsen")


if __name__ == "__main__":
    unittest.main()
