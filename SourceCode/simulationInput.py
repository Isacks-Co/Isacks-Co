from dataclasses import dataclass
from ase import Atoms

class SimulationSettings:
    """
    Abstract class for simulation settings.
    Shouldn't be used by itself
    """
    def __init__(self, timestep: float, total_time: float, potential: str):
        self.timestep = timestep
        self.total_time = total_time
        self.potential = potential


class NVESettings(SimulationSettings):
    @property
    def ensemble(self):
        return "NVE"

    def __str__(self):
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Potential : {self.potential} \n"
            f"Timestep: {self.timestep} fs  \n"
            f"Total time : {self.total_time} fs \n"
        )


class NVTSettings(SimulationSettings):
    def __init__(self, timestep: float, total_time: float, potential: str, temperature: float):
        super().__init__(timestep, total_time, potential)
        self.temperature = temperature

    @property
    def ensemble(self):
        return "NVT"

    def __str__(self):
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Potential : {self.potential} \n"
            f"Temperature: {self.temperature} K\n"
            f"Timestep: {self.timestep} fs \n"
            f"Total time : {self.total_time} fs \n"
        )


class NPTSettings(SimulationSettings):
    def __init__(self, timestep: float, total_time: float, potential: str, temperature: float, pressure: float):
        super().__init__(timestep, total_time, potential)
        self.temperature = temperature
        self.pressure = pressure

    @property
    def ensemble(self):
        return "NPT"

    def __str__(self):
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Potential : {self.potential} \n"
            f"Temperature: {self.temperature} K\n"
            f"Pressure: {self.pressure} Gpa \n"
            f"Timestep: {self.timestep} fs \n"
            f"Total time : {self.total_time} fs \n"
        )


class SimulationInput:
    def __init__(self, atoms, settings: SimulationSettings):
        self.atoms = atoms
        self.settings = settings

    def __str__(self):
        return (
            f"Atomic structure: {self.atoms} \n"
            f"{self.settings}"
        )



if __name__ == "__main__":
    settings = NPTSettings(timestep=5,total_time=250,potential="EMT",temperature=200,pressure=20)
    simInput = SimulationInput("a",settings)
    


