import sys

sys.path.append("..")
from Tests.TestBase import TestBase
from SourceCode.PreProcessing import PreProcessing
from ase.visualize import view

class TestPreProcessing(TestBase):
    """Tests for class PreProcessing"""
    def setUp(self):
        super().setUp()
        self.Pre = PreProcessing("TestSettings/solidSettings.json", "TestAtomicStructure/Cu_fcc.vasp", None)
    
    def testSettings(self):
        # Make sure ValueError is raised for bad settings.
        self.Pre.settings["Temperature"] = 4000

        with self.assertRaises(ValueError):
            self.Pre.sanityCheckSettings()

    def testAtomicStructure(self):
        """Test that supercells is working properly."""
        assert len(self.Pre.atoms) == 125

    def testCreateMD(self):
        MD = self.Pre.createMD()
        assert str(MD.integrator).startswith("functools.partial(<function Langevin")

    def testTerminalRead(self):
        self.Pre.readTerminalInput(["-T", 50])
        assert self.Pre.settings["Temperature"] == 50