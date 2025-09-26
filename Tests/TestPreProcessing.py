import sys

sys.path.append("..")
from Tests.TestBase import TestBase
from PreProcessing import PreProcessing
from ase.visualize import view

class TestPreProcessing(TestBase):

    def testSettings(self):
        PP = PreProcessing("testsettings.json", "TESTPOSCAR", None)
        # Make sure ValueError is raised for bad settings.
        PP.settings["Temperature"] = 4000

        with self.assertRaises(ValueError):
            PP.sanityCheckSettings()

    def testAtomicStructure(self):
        """Test that supercells is working properly."""
        PP = PreProcessing("testsettings.json", "TESTPOSCAR", None)
        assert len(PP.atoms) == 125

    def testCreateMD(self):
        PP = PreProcessing("testsettings.json", "TESTPOSCAR", None)
        MD = PP.createMD()
        assert str(MD.integrator).startswith("functools.partial(<function VelocityVerlet")