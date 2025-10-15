import sys

sys.path.append("..")
import logging
import pytest
from Tests.TestBase import TestBase
from MolecularDynamics import main as MCMain

logger = logging.getLogger(__name__)


class TestMolecularDynamics(TestBase):
    """Run premade cases that are used in TestPostProcessing later. Order is important!"""
    def setUp(self):
        super().setUp()
        sys.argv = [sys.argv[0]]

    def testMeltedCu(self):
        MCMain(settings="TestSettings/meltedSettings.json", poscar="TestAtomicStructure/Cu_fcc.vasp")

    def testSolidCu(self):
        MCMain(settings="TestSettings/solidSettings.json", poscar="TestAtomicStructure/Cu_fcc.vasp")

    def testNearZeroCu(self):
        MCMain(settings="TestSettings/nearZeroSettings.json", poscar="TestAtomicStructure/Cu_fcc.vasp")

    def testBccCu(self):
        MCMain(settings="TestSettings/chromiumSettings.json", poscar="TestAtomicStructure/Cr_bcc.vasp")

    def testScCu(self):
        MCMain(settings="TestSettings/poloniumSettings.json", poscar="TestAtomicStructure/Po_sc.vasp")