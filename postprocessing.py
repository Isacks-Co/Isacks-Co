from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.io import read
from asap3.visualize.primiplotter import *
from numpy import *

class PostProcessing:
    """
        Uses the log files generates from the Molecular Dynamics class to vizualize and analyze the data
        input trajectory file as traj_file 
    """
    def __init__(self, traj_file):
        self.read_traj_file = traj_file
        
    def vizualize(self):
        for atoms in self.read_traj_file:
            plotter = PrimiPlotter(atoms)
            plotter.set_output(X11Window())   # Plot in X11 window on the user's screen
        #view(read(self.read_traj_file), ':')
    