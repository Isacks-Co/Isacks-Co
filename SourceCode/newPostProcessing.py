import logging
from simulationInput import SimulationSettings
from ase.io.trajectory import Trajectory
from quantityCalculator import QuantityCalculator
from Utils import plot
from scipy.interpolate import UnivariateSpline
import numpy as np


class PostProcessing():
    def __init__(self, settings : SimulationSettings):
        self.settings = settings
        self.traj = Trajectory(f"{self.settings.output_file}.traj")
        
        equil_traj = Trajectory("../Outputs/equil_output_file.traj")
        self.qce = QuantityCalculator(self.settings, equil_traj)
        self.qc = QuantityCalculator(self.settings, self.traj)

        """
        msds =[]
        e_pot = []
        e_kin = []
        x= []
        for frame in range(len(equil_traj)):
            x.append(frame)
            msd = self.qc.computeMSD(frame)
            msds.append(msd)
            e_pot.append(equil_traj[frame].get_potential_energy())
            #e_kin.append(frame.get_kinetic_energy())
        x_new = np.linspace(0,max(x),500)
        spline_epot = UnivariateSpline(x,e_pot,s = 0.5)(x_new)

        print(spline_epot)
        plot(spline_epot)
        #plot(e_pot)
        #plot(e_kin)
        """
