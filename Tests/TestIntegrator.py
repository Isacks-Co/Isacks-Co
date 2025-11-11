import sys
import numpy as np
sys.path.append("..")
from Tests.TestBase import TestBase
from SourceCode.ASEWrappers.integrator import Integrator,VelocityVerletIntegrator,LangevinIntegrator,IsotropicMTKNPTIntegrator,BerendsenNPTIntegrator

import unittest

import logging

logger = logging.getLogger(__name__)

class TestIntegrator(TestBase):
    """Tests for the class Integrator class. """
    def generalTest(self):
        integrator = Integrator(1)
        assert len(integrator.attachments) == 0
        integrator.attach(1,1)
        assert len(integrator.attachments) == 1

    def testVelocityVerlet(self):
        integrator = VelocityVerletIntegrator(1)
        assert integrator.ensemble == "NVE"
        
        assert  str(integrator.integrator_partial).startswith("functools.partial(<function VelocityVerlet")
        
    def testLangevin(self):
        integrator = LangevinIntegrator(1,300,0.01)

        assert integrator.ensemble == "NVT"
        print(integrator.integrator_partial)
        assert  str(integrator.integrator_partial).startswith("functools.partial(<function Langevin")
        
    def testIsotropicMTK(self):
        integrator = IsotropicMTKNPTIntegrator(1,300,0,1,1)
        assert integrator.ensemble == "NPT"
        print(integrator.integrator_partial)
        assert  str(integrator.integrator_partial).startswith("functools.partial(<class 'asap3.md.nose_hoover_chain.IsotropicMTKNPT'>")
        
    def testNPTBerendsen(self):
        integrator = BerendsenNPTIntegrator(1,300,0,1e-7)
        assert integrator.ensemble == "NPT"
        print(integrator.integrator_partial)
        assert  str(integrator.integrator_partial).startswith("functools.partial(<class 'asap3.md.nptberendsen")
if __name__ == "__main__":
    unittest.main()