import sys

sys.path.append("../SourceCode")
import logging
import pytest
from TestBase import TestBase
from MolecularDynamics import main as MolecularDynamics

logger = logging.getLogger(__name__)


class TestMolecularDynamics(TestBase):
    """Run premade cases that are used in TestPostProcessing later. Order is important!"""
    def setUp(self):
        super().setUp()

    def testMeltedCu(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/meltedSettings.json"]
        MolecularDynamics()

    def testSolidCu(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/solidSettings.json"]
        MolecularDynamics()

    def testNearZeroCu(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/nearZeroSettings.json"]
        MolecularDynamics()
    @pytest.mark.skip("Lennard Jones doesnt work")
    def testBccCu(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cr_bcc.vasp", "TestSettings/chromiumSettings.json"]
        MolecularDynamics()

    @pytest.mark.skip("NPT currently not working. Half implemented, something went wrong during merge")
    def testNPTCopper(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/NPTCopperSettings.json"]
        MolecularDynamics()