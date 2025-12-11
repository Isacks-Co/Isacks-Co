# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import logging
import pytest
import sys
sys.path.append("../SourceCode")
from MolecularDynamics import main as MolecularDynamics
from TestBase import TestBase
import subprocess
import pathlib
import unittest



logger = logging.getLogger(__name__)


class TestMolecularDynamics(TestBase):
    """Run some simulations. Was previously used in TestQuantityCalculator but not anymore.
    May be used in TestPostProcessing later on or removed
    NOTE: Currently out of the CI action on github"""
    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)
        self.BASE_DIR = pathlib.Path(__file__).parent
        self.SIMULATION_DIR = self.BASE_DIR / "TestSimulation"
        self.script_name = "runMD.sh"
    def setUp(self):
        super().setUp()

    def testFullProgram(self):
        try:
            command = ["./" + self.script_name, "settings.json"]
            subprocess.run(command,
            cwd=self.SIMULATION_DIR,  # Execute the command as if launched from this directory
            capture_output=True,
            text=True,
            check=True)
        
        except subprocess.CalledProcessError as e:
            pytest.fail(
            f"Simulation script failed with exit code {e.returncode}.\n"
            f"STDOUT:\n{e.stdout}\n"
            f"STDERR:\n{e.stderr}"
        )
        
    # If the code reaches here, the script exited with 0, and the test passes.
    assert True
if __name__ == "__main__":
    unittest.main()


