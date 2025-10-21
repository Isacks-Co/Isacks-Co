import logging
from asap3 import LennardJones as asap_LJ
from asap3 import EMT
from mace.calculators import MACECalculator
import torch
import warnings

log = logging.getLogger(__name__)


class Potential:

    def __init__(self):
        self.pot_str = None

    def getPotentialCalculator(self):
        return None

    def __str__(self):
        return self.pot_str


class LennardJonesPotential(Potential):

    def __init__(self, sigma: float, eps: float, rc: float = None, ro: float = None):
        super().__init__()
        self.pot_str = "Lennard Jones"
        self.sigma = sigma
        self.eps = eps
        self.rc = rc if rc is not None else 2.5 * sigma
        self.ro = ro if ro is not None else 0.9 * rc

    def getPotentialCalclulator(self):
        return asap_LJ(sigma=self.sigma, eps=self.eps, rc=self.rc, ro=self.ro)


class EMTPotential(Potential):

    def __init__(self):
        super().__init__()
        self.pot_str = "EMT"

    def getPotentialCalculator(self):
        return EMT()


class MACEPotential(Potential):

    def __init__(self, model_path):
        super().__init__()
        self.model_path = model_path

    def getPotentialCalculator(self):
        device = "cuda" if torch.cuda.is_available() else "cpu"
        return MACECalculator(model_paths=self.model_path , device =device, head="default")




