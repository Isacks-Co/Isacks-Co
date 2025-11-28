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
