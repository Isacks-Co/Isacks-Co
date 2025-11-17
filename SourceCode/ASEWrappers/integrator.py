from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from asap3.md.nose_hoover_chain import IsotropicMTKNPT
from asap3.md.nptberendsen import NPTBerendsen
from ase.units import fs,GPa

from functools import partial
from copy import deepcopy

from .contextManager import suppress_cpp_output

"""
Wrapper of the ASE integrators. 
The purpose is to make the code more readable 
by wrapping up some of the complexities ASEs 
interface results in. 
"""


class Integrator:
    def __init__(self,timestep):
        self.timestep = timestep
        self.atoms = None
        self.integrator_partial = None
        self.attachments = []
    
    @property
    def ensemble(self):
        return None

    def attach(self,attachment,interval):
        self.attachments.append( (attachment,interval))

    def _addAttachments(self,integrator):
        for func,interval in self.attachments:

            integrator.attach(func,interval = interval)
        
    def run(self,atomic_structure,steps):
        with suppress_cpp_output():
            integrator_func = deepcopy(self.integrator_partial)(atomic_structure.getAtoms())

            self._addAttachments(integrator_func)

            integrator_func.run(steps)
        self.clearData()
       

    def clearData(self):
        self.atoms = None
        self.attachments = []
    def __str__(self): # TODO Expand on this
        return self.ensemble
    
class VelocityVerletIntegrator(Integrator):
    def __init__(self,timestep):
        super().__init__(timestep=timestep)
        self.integrator_partial = partial(VelocityVerlet, timestep=self.timestep *fs)
    
    @property
    def ensemble(self):
        return "NVE"
    
class LangevinIntegrator(Integrator):
    def __init__(self,timestep,temperature_K,friction):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.friction = friction 
        self.integrator_partial = partial(Langevin, timestep=self.timestep*fs, temperature_K=self.temperature_K,
                                     friction=self.friction/fs)
    @property
    def ensemble(self):
        return "NVT"
    
class IsotropicMTKNPTIntegrator(Integrator):
    def __init__(self,timestep,temperature_K,pressure,tdamp,pdamp):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.tdamp = tdamp
        self.pdamp = pdamp
        self.integrator_partial = partial(IsotropicMTKNPT, timestep=self.timestep*fs, temperature_K=self.temperature_K,
                                     pressure_au=self.pressure * GPa, tdamp = self.tdamp,pdamp = self.pdamp)
        
    @property
    def ensemble(self):
        return "NPT"
    

class BerendsenNPTIntegrator(Integrator):
    def __init__(self,timestep,temperature_K,pressure,compressibility):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.compressibility = compressibility
        self.integrator_partial = partial(NPTBerendsen, timestep=self.timestep*fs, temperature_K=self.temperature_K,
                                     pressure_au=self.pressure * GPa,compressibility = self.compressibility / GPa)
        
    @property
    def ensemble(self):
        return "NPT"