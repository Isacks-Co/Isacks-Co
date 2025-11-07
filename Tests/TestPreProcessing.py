import sys

sys.path.append("../SourceCode")
from TestBase import TestBase
from PreProcessing import PreProcessing
from ase.visualize import view

class TestPreProcessing(TestBase):
    """Tests for class PreProcessing"""
    def setUp(self):
        super().setUp()
        sys.argv = [sys.argv[0], "TestAtomicStructure/Cu_fcc.vasp", "TestSettings/meltedSettings.json"]
        self.Pre = PreProcessing(sys.argv)
    
    def testSettings(self):
        # Make sure ValueError is raised for bad settings.
        self.Pre.settings["Temperature"] = -1

        with self.assertRaises(ValueError):
            self.Pre.sanityCheckSettings()

    def testAtomicStructure(self):
        """Test that supercells is working properly."""
        assert len(self.Pre.atoms) == 125


    def testTerminalRead(self):
        sys.argv = sys.argv + ["-T", "50"]
        self.Pre = PreProcessing(sys.argv)
        assert self.Pre.createSettings().temperature == 50