from asap3 import LennardJones as asap_LJ
from asap3 import EMT



class Potential:
    
    def __init__(self):
        self.pot_str = None

    def getASEPotentialCalculator(self):
        return None

    def __str__(self):
        return self.pot_str


class LennardJonesPotential(Potential):

    def __init__(self, sigma: float, eps: float, rc: float = None, ro: float = None):
        super().__init__()
        self.pot_str = "Lennard Jones"
        self.sigma = sigma
        self.eps = eps
        self.rc = rc if rc is not None else 2.5*sigma
        self.ro = ro if ro is not None else 0.9*rc

    def getASEPotentialCalculator(self):
        return asap_LJ(sigma = self.sigma,eps = self.eps,rc = self.rc, ro = self.ro)

class EMTPotential(Potential):

    def __init__(self):
        super().__init__()
        self.pot_str = "EMT"         

    def getASEPotentialCalculator(self):
        return EMT()
