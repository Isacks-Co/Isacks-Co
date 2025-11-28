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
import sys
sys.path.append("../SourceCode")
from TestBase import TestBase
from Utils.LJRegistry import LJParams

logger = logging.getLogger(__name__)


class TestLJParams(TestBase):

    def testNonEmptyArgs(self):
        logger.info("Testing setting dictionary with empty sigma and epsilon parameters")
        known_material = LJParams(material="Ar")
        self.assertIsInstance(known_material, dict)
        self.assertGreater(known_material["epsilon_eV"], 0.0)
        self.assertGreater(known_material["sigma_A"], 0.0)

    def testRadiusArgs(self):
        known_material = LJParams(material="Ar")
        self.assertLess(known_material["ro_A"], known_material["rc_A"])

    def testCaseInsensitive(self):
        logger.info("Testing case insesitivity")
        self.assertEqual(LJParams("ar"), LJParams("Ar"))
        self.assertEqual(LJParams("AR"), LJParams("aR"))
        self.assertEqual(LJParams("ar"), LJParams("AR"))

    def testOverride(self):
        logger.info("Checking that the override values are used instead of the ones from database")
        known = LJParams(material="Ar")
        override = LJParams(material="Ar", epsilon_eV=0.02, sigma_A=3.5)
        self.assertNotEqual(known["epsilon_eV"], override["epsilon_eV"])
        self.assertNotEqual(known["sigma_A"], override["sigma_A"])


if __name__ == "__main__":
    import unittest

    unittest.main(verbosity=2)
