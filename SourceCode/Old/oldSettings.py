from SourceCode.simulationInput import SimulationSettings
#TODO Remove.  Ensemble is directly related to the choice of integrator


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
        parent_string = super().__str__()
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Initial temperature: {self.initial_temperature} K \n"
            f"{parent_string}\n"
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
        parent_string = super().__str__()
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Temperature: {self.temperature} K\n"
            f"Friction: {self.friction}\n"
            f"{parent_string} \n"
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
        float tdamp: Time constant for temperature damping
        float pdamp: Time constant for pressure damping

    """
    def __init__(self,timestep: float, num_steps: int, potential: str, interval: int, temperature: float,pressure: float, tdamp: float, pdamp: float, supercells: list, output_file: str = "output.traj"):
        super().__init__(timestep, num_steps, potential, output_file, supercells, interval)
        self.temperature = temperature
        self.pressure = pressure
        self.tdamp = tdamp
        self.pdamp = pdamp

    @property
    def ensemble(self):
        return "NPT"

    def __str__(self):
        parent_string = super().__str__()
        return (
            f"Ensemble: {self.ensemble} \n"
            f"Temperature: {self.temperature} K\n"
            f"Pressure: {self.pressure} Pa \n"
            f"{parent_string} \n"

        )


