import sys

sys.path.append("..")
from Tests.TestBase import TestBase
from PreProcessing import PreProcessing
from ase.visualize import view

class TestPreProcessing(TestBase):

    def testSettings(self):
        PP = PreProcessing("Tests/testsettings.json", "Tests/TESTPOSCAR", None)
        PP.sanityCheckSettings()

    def testAtomicStructure(self):
        PP = PreProcessing("Tests/testsettings.json", "Tests/TESTPOSCAR", None)
        assert len(PP.atoms) == 125

    def testCreateMD(self):
        PP = PreProcessing("Tests/testsettings.json", "Tests/TESTPOSCAR", None)
        MD = PP.createMD()
        assert str(MD.integrator).startswith("functools.partial(<function VelocityVerlet")