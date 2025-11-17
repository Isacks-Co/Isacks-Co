import logging
from simulationInput import SimulationSettings
from ase.io.trajectory import Trajectory
from quantityCalculator import QuantityCalculator


class PostProcessing():
    def __init__(self, settings : SimulationSettings):
        self.settings = settings
        self.traj = Trajectory(f"{self.settings.output_file}.traj")
        self.qc = QuantityCalculator(self.settings, self.traj)
        self.qc.getQuantities()

