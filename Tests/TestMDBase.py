import sys

sys.path.append("..")
from Tests.TestBase import TestBase
from SourceCode.MDBase import MDBase
from ase.lattice.cubic import FaceCenteredCubic
import unittest
from asap3 import EMT as asap_EMT
import logging

logger = logging.getLogger(__name__)


class TestMDBase(TestBase):
    """Tests for the class MDBase"""

    def testNVEinit(self):
        """Check that the correct function is attached for all the ensembles"""
        MD = MDBase.initNVE(temperature=293, pot_str="EMT", timestep=0.5, steps=500, interval=5, output_file="data",
                            equilibrium_steps=5000)
        logger.info("Checking that NVE has VelocityVerlet")
        assert str(MD.integrator).startswith("functools.partial(<function VelocityVerlet")

    def testNVTinit(self):
        MD = MDBase.initNVT(temperature=293, friction=0.01, pot_str="EMT", timestep=0.5, steps=500, interval=5,
                            output_file="data", equilibrium_steps=5000)
        assert str(MD.integrator).startswith("functools.partial(<function Langevin")

    def testNPTinit(self):
        MD = MDBase.initNPT(temperature=293, timestep=0.5, steps=500, interval=5, pressure_Pa=10e+6,
                            compressibility=10e-11, pot_str="EMT", output_file="data", equilibrium_steps=5000)

        assert str(MD.integrator).startswith("functools.partial(<class 'asap3.md.nptberendsen")

    def testEquilibrium(self):
        atoms = FaceCenteredCubic(size=(5, 5, 5), symbol="Cu", pbc=True)
        atoms.calc = asap_EMT()
        MD = MDBase()
        MD.equilibriumRun(atoms)
        # TODO implement good test to see that we reach equil.


if __name__ == "__main__":
    unittest.main()
