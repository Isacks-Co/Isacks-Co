import sys

sys.path.append("../SourceCode")
<<<<<<< HEAD
from TestBase import TestBase
from MDBase import MDBase
from PreProcessing import PreProcessing
=======
from Tests.TestBase import TestBase
from SourceCode.Old.MDBase import MDBase
>>>>>>> aa47b38 (Changed some of the workflow. Now runnig using scripts in Simulation folder. Currently supporting setup, full run and clearing. Will add specific ones for each part of the full run like equilRun and Postprocessing)
from ase.lattice.cubic import FaceCenteredCubic
import pytest
import sys
from asap3 import EMT as asap_EMT
import logging

logger = logging.getLogger(__name__)


class TestMDBase(TestBase):
    """Tests for the class MDBase"""

    def setUp(self):
        super().setUp()
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/solidSettings.json"]
        pre_processing = PreProcessing(sys.argv)
        self.MD = MDBase(pre_processing.createSettings())

    def testNVEinit(self):
        """Check that the correct function is attached for all the ensembles"""
        logger.info("Checking that NVE has VelocityVerlet")
        assert str(self.MD.getIntegrator("NVE")).startswith("functools.partial(<function VelocityVerlet")

    def testNVTinit(self):
        assert str(self.MD.getIntegrator("NVT")).startswith("functools.partial(<function Langevin")

    @pytest.mark.skip("NPT not working currently")
    def testNPTinit(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/NPTCopperSettings.json"]
        pre_processing = PreProcessing(sys.argv)
        MDNPT = MDBase(pre_processing.createSettings())

        assert str(MDNPT.getIntegrator("NPT")).startswith("functools.partial(<class 'asap3.md.nptberendsen")

