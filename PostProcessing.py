from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.io import read

class PostProcessing:
    """
        Uses the log files generates from the Molecular Dynamics class to vizualize and analyze the data
        input trajectory file as traj_file 
    """
    def __init__(self, traj_file):
        self.read_traj_file = traj_file
        
    def vizualize(self):
        view(read(self.read_traj_file), ':')