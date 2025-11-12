from asap3 import LennardJones as asap_LJ
from asap3 import EMT

#from mace.calculators import MACECalculator
#import torch

class Potential:
    
    def __init__(self):
        self.pot_str = None

    def getASEPotentialCalculator(self):
        return None

    def __str__(self):
        return self.pot_str


class LennardJonesPotential(Potential):

    def __init__(self,atomic_numbers, sigmas, epsilons, rc: float = None, ro: float = None):
        super().__init__()
        self.pot_str = "Lennard Jones"
        self.atomic_nums = atomic_numbers
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.rc = rc if rc is not None else 2.5*max(self.sigmas)
        #self.ro = ro if ro is not None else 0.9*self.rc

    def getASEPotentialCalculator(self):
        return asap_LJ(self.atomic_nums, sigma = self.sigmas,epsilon = self.epsilons,rCut = self.rc)

class EMTPotential(Potential):

    def __init__(self):
        super().__init__()
        self.pot_str = "EMT"         

    def getASEPotentialCalculator(self):
        return EMT()
    

"""
class MACEPotential(Potential):
    def __init__(self,model_path = "MACEModels/mace-mpa-0-medium.model"):
        super().__init__()
        self.pot_str = "MACE"  
        self.model_path = model_path

    def getASEPotentialCalculator(self):
        device = "cuda" if torch.cuda.is_available() else "cpu"
        return MACECalculator(model_paths= self.model_path, device=device, default_dtype="float64",head="default")
"""