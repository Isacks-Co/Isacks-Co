import sys

sys.path.append("../SourceCode")
import logging
import pytest
from TestBase import TestBase
from MolecularDynamics import main as MolecularDynamics

logger = logging.getLogger(__name__)


class TestMolecularDynamics(TestBase):
    """Run some simulations. Was previously used in TestQuantityCalculator but not anymore.
    May be used in TestPostProcessing later on or removed
    NOTE: Currently out of the CI action on github"""

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

    def testBccCu(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cr_bcc.vasp", "TestSettings/chromiumSettings.json"]
        MolecularDynamics()

    def testNPTCopper(self):
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/NPTCopperSettings.json"]
        MolecularDynamics()
