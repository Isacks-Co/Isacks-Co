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


import os
import torch
from asap3 import EMT
from asap3 import LennardJones as asap_LJ
from mace.calculators import MACECalculator


class Potential:

    def __init__(self):
        self.pot_str = None

    def getASEPotentialCalculator(self):
        return None

    def __str__(self):
        return self.pot_str


class LennardJonesPotential(Potential):

    def __init__(self, atomic_numbers, sigmas, epsilons, rc=None):
        super().__init__()
        self.pot_str = "Lennard Jones"

        self.atomic_nums = atomic_numbers
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.rc = rc if rc is not None else 2.5 * max(self.sigmas)

        # self.ro = ro if ro is not None else 0.9*self.rc

    def getASEPotentialCalculator(self):
        return asap_LJ(self.atomic_nums, sigma=self.sigmas, epsilon=self.epsilons, rCut=self.rc)


class EMTPotential(Potential):

    def __init__(self):
        super().__init__()
        self.pot_str = "EMT"

    def getASEPotentialCalculator(self):
        return EMT()


class MACEPotential(Potential):
    def __init__(self, model_path=os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))) + "/MACEModels/mace-mpa-0-medium.model"):
        super().__init__()
        self.pot_str = "MACE"
        self.model_path = model_path

    def getASEPotentialCalculator(self):
        import os, warnings
        os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
        warnings.filterwarnings(
            "ignore",
            message="Environment variable TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD detected, since the`weights_only` argument was not explicitly passed to `torch.load`, forcing weights_only=False.",
        )
        device = "cuda" if torch.cuda.is_available() else "cpu"
        return MACECalculator(model_paths=self.model_path, device=device, default_dtype="float64", head="default")
