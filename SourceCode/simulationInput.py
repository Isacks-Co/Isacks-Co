from dataclasses import dataclass
from ase import Atoms


@dataclass
class SimulationSettings:
    """
    Abstract class for simulation settings.
    Shouldnt be used by itself
    """
    timestep: float
    total_time: float
    potential: str

@dataclass
class NVESettings(SimulationSettings):

    @property
    def ensemble(self): 
        return "NVE"
    
    def __str__(self):
        return  (
                f"Ensemble: {self.ensemble} \n"
                f"Potential : {self.potential} \n"
                f"Timestep: {self.timestep} fs  \n"
                f"Total time : {self.total_time} fs \n"
                )

@dataclass
class NVTSettings(SimulationSettings):
    temperature: float
    
    @property
    def ensemble(self): 
        return "NVT"

    def __str__(self):
        return  (
                f"Ensemble: {self.ensemble} \n"
                f"Potential : {self.potential} \n"
                f"Temperature: {self.temperature } K\n"
                f"Timestep: {self.timestep} fs \n"
                f"Total time : {self.total_time} fs \n"
                )
     

@dataclass
class NPTSettings(SimulationSettings):
    temperature: float
    pressure: float
    
    @property
    def ensemble(self): 
        return "NPT"
    
    def __str__(self):
        return  (
                f"Ensemble: {self.ensemble} \n"
                f"Potential : {self.potential} \n"
                f"Temperature: {self.temperature } K\n"
                f"Pressure: {self.pressure } Gpa \n"
                f"Timestep: {self.timestep} fs \n"
                f"Total time : {self.total_time} fs \n"

                )
    

@dataclass
class SimulationInput:
    atoms: Atoms 
    settings: SimulationSettings

    def __str__(self) -> str:
        return  (
                f"Atomic structure: {self.atoms} \n"
                f"{self.settings}"
                )


if __name__ == "__main__":
    settings = NPTSettings(timestep=5,total_time=250,potential="EMT",temperature=200,pressure=20)
    simInput = SimulationInput("a",settings)
    


