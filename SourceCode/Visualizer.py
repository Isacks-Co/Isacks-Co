from ase.visualize import view
from ase.io.trajectory import Trajectory

class Visualizer():
    def __init__(self, traj : Trajectory):
        self.traj = traj

    def visualize(self):
        view(self.traj)