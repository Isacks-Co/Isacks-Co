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
import unittest
sys.path.append("../SourceCode")
from ASEWrappers import LennardJonesPotential, EMTPotential
from TestBase import TestBase
from asap3.Internal.BuiltinPotentials import EMT
from asap3.Internal.BuiltinPotentials import LennardJones


logger = logging.getLogger(__name__)


class TestPotential(TestBase):
    """Tests for the class Potential class. """

    def test_EMT(self):
        pot = EMTPotential()
        assert str(pot) == "EMT"

        assert isinstance(pot.getASEPotentialCalculator(), EMT)

    def test_LJ(self):
        pot = LennardJonesPotential(atomic_numbers=1, sigmas=[0.5], epsilons=[0.5])

        assert str(pot) == "Lennard Jones"

        assert isinstance(pot.getASEPotentialCalculator(), LennardJones)


if __name__ == "__main__":
    unittest.main()
