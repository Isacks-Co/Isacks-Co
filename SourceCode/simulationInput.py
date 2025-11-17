from ASEWrappers import Potential,Integrator

class SimulationSettings:
    """
    args:
        
        (int) num_step: Number of steps in the simulation
        (Potential) potential : Potenial object defining the interaction potential
        (Integrator) integrator : Integrator object defining how the system should be integrated through time
        
    """
    def __init__(self, num_steps: int, potential: Potential, integrator: Integrator ):
        
        self.num_steps = num_steps
        self.integrator = integrator
        self.potential = potential

    def __str__(self):
        return (
            f"Number of steps: {self.num_steps} fs \n"
            f"Potential: {self.potential} \n"
            f"Integrator: {self.integrator}\n"
        )
        
