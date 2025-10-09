from ase import Atoms

class SimulationSettings:
    """
    Abstract class for simulation settings.
    Shouldn't be used by itself
    args:
        float timestep: size of the timesteps in fs
        int num_step: Number of steps in the simulation
        int sample_interval: Size of the interval we sample simulation for data such as energy, temperature etc. Given as a number of timesteps
        str potential: Potential to use for the simulations such as Lennard Jones, EMT ...
        str output_file: Path to write output --> output_file.traj
        
    """
    def __init__(self, timestep: float, num_steps: int, potential: str, output_file: str , supercells: list, interval: int = 5 ):
        self.timestep = timestep
        self.num_steps = num_steps
        self.sample_interval = interval
        self.potential = potential
        self.output_file = output_file
        self.supercells = supercells
        

class NVESettings(SimulationSettings):
    """
    SimulationSettings object holding the specific parameters for NVE simulations. Is mainly used
    as a container for easy acces to params.


    args:
        float timestep: size of the timesteps in fs
        int num_step: Number of steps in the simulation
        int sample_interval: Size of the interval we sample simulation for data such as energy, temperature etc. Given as a number of timesteps
        str potential: Potential to use for the simulations such as Lennard Jones, EMT ...
        str output_file: Path to write output --> output_file.traj

        float init_temp: Temperature to initialize the structure at in Kelvin. Basically directly related
                        to the kinetic energy. 
    """

    def __init__(self,timestep: float, num_steps: int, potential: str, interval: int, init_temp: float, supercells: list, output_file: str = "output.traj"):
        super().__init__(timestep, num_steps, potential, output_file, supercells, interval)
        self.initial_temperature = init_temp

    @property
    def ensemble(self):
        return "NVE"

    def __str__(self):
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Potential : {self.potential} \n"
            f"Timestep: {self.timestep} fs  \n"
            f"Number of steps: {self.num_steps} \n"
        )


class NVTSettings(SimulationSettings):

    """
    SimulationSettings object holding the specific parameters for NVT simulations. Is mainly used
    as a container for easy acces to params.


    args:
        float timestep: size of the timesteps in fs
        int num_step: Number of steps in the simulation
        int sample_interval: Size of the interval we sample simulation for data such as energy, temperature etc. Given as a number of timesteps
        str potential: Potential to use for the simulations such as Lennard Jones, EMT ...
        str output_file: Path to write output --> output_file.traj

        float temperature:  Temperature to simulate the structure at in Kelvin. Note difference between this and init_temperature.
                            Here we want to keep this temperature constant. Where in NVE we let it change and only needed a starting point.
        float friction: Friction for simulatio in fs^-1. Controlls how quickly temperature changes. Large friction -> less change
    """
    def __init__(self,timestep: float, num_steps: int, potential: str, interval: int, temperature: float, friction: float, supercells: list,  output_file: str = "output.traj"):
        super().__init__(timestep, num_steps, potential, output_file, supercells, interval )
        self.temperature = temperature
        self.friction = friction 

    @property
    def ensemble(self):
        return "NVT"

    def __str__(self):
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Potential : {self.potential} \n"
            f"Temperature: {self.temperature} K\n"
            f"Timestep: {self.timestep} fs \n"
            f"Number of steps : {self.num_steps} \n"
        )


class NPTSettings(SimulationSettings):
    """
    SimulationSettings object holding the specific parameters for NVT simulations. Is mainly used
    as a container for easy acces to params.


    args:
        float timestep: size of the timesteps in fs
        int num_step: Number of steps in the simulation
        int sample_interval: Size of the interval we sample simulation for data such as energy, temperature etc. Given as a number of timesteps
        str potential: Potential to use for the simulations such as Lennard Jones, EMT ...
        str output_file: Path to write output --> output_file.traj

        float temperature:  Temperature to simulate the structure at in Kelvin. Note difference between this and init_temperature.
                            Here we want to keep this temperature constant. Where in NVE we let it change and only needed a starting point.
        float pressure: Pressure to simulate the structure at in Pa 
        float compressibility: Compressibility for the NPT simulation in units Pa^-1. Large comp -> box volume fluctuates easily

    """
    def __init__(self,timestep: float, num_steps: int, potential: str, interval: int, temperature: float,pressure: float, compressibility: float, supercells: list, output_file: str = "output.traj"):
        super().__init__(timestep, num_steps, potential, output_file, supercells, interval)
        self.temperature = temperature
        self.pressure = pressure
        self.compressibility = compressibility

    @property
    def ensemble(self):
        return "NPT"

    def __str__(self):
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Potential : {self.potential} \n"
            f"Temperature: {self.temperature} K\n"
            f"Pressure: {self.pressure} Pa \n"
            f"Timestep: {self.timestep} fs \n"
            f"Number of steps : {self.num_steps} \n"
        )




if __name__ == "__main__":
    settings = NPTSettings(timestep=5,num_steps=250,potential="EMT",temperature=200,pressure=20,interval=2,compressibility=1)
   


