from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from asap3.md.nose_hoover_chain import IsotropicMTKNPT
from asap3.md.nptberendsen import NPTBerendsen
from ase.units import fs,GPa
from functools import partial

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

    def _addAttachments(self,integrator,atomic_structure):
        for func,interval in self.attachments:
            integrator.attach(lambda: func(atomic_structure),interval = interval)

    def run(self,atomic_structure,steps):
        integrator_func = self.integrator_partial(atomic_structure.get_atoms())
        integrator_func = self._addAttachments(integrator_func,atomic_structure)
        self.integrator_func.run(steps)

    def __str__(self):
        return self.ensemble
    
class VelocityVerletIntegrator(Integrator):
    def __init__(self,timestep):
        super().__init__(timestep=timestep)
        self.partial = partial(VelocityVerlet, timestep=self.timestep *fs)
    
    @property
    def ensemble(self):
        return "NVE"
    
class LangevinIntegrator(Integrator):
    def __init__(self,timestep,temperature_K,friction):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.friction = friction 
        self.integrator_func = partial(Langevin, timestep=self.timestep*fs, temperature_K=self.temperature_K,
                                     friction=self.friction/fs)
    @property
    def ensemble(self):
        return "NVT"
    
class IsoTropicMTKNPTIntegrator(Integrator):
    def __init__(self,timestep,temperature_K,pressure,tdamp,pdamp):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.tdamp = tdamp
        self.pdamp = pdamp
        self.integrator_func = partial(IsotropicMTKNPT, timestep=self.timestep*fs, temperature_K=self.temperature_K,
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
        self.integrator_func = partial(NPTBerendsen, timestep=self.timestep*fs, temperature_K=self.temperature_K,
                                     pressure_au=self.pressure * GPa,compressibility = self.compressibility / GPa)
        
    @property
    def ensemble(self):
        return "NPT"