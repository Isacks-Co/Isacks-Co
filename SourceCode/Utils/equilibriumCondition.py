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
import numpy as np

log = logging.getLogger(__name__)  # Nothing will prob log here?


class EquilibriumCondition:

    @staticmethod
    def checkInternalPressureStable(list_of_internal_pressure, tol=1e-3, ):
        """
        Used to check if the difference in internal pressure is below threshold.

        Returns:
            bool: True if the difference is below threshold; False otherwise.
        """

        mean_pressure = np.mean(list_of_internal_pressure)
        print(mean_pressure)
        return  (np.isclose(mean_pressure, 0, atol=tol))

    @staticmethod
    def checkStable(list_of_values, tol=1e-5):
        """
        Used to check if energy oscillates around certain value, then we probably found the equilibrium,
        meant to input either a list of MSD values or energy, which unlike internal pressure don't have to oscillate
        around 0.

        return:
            true if the difference is below threshold
            false if the difference is above threshold
        """

        delta = np.abs(np.max(list_of_values) - np.min(list_of_values))

        return delta <= tol
