from ase.io.trajectory import Trajectory
from simulationInput import SimulationSettings
import numpy as np
import logging

logger = logging.getLogger(__name__)

class QuantityCalculator:

    def __init__(self,settings: SimulationSettings, traj : Trajectory):
        self.traj = traj
        self.settings = settings


    def getQuantities(self):
        #  Compute all general quantities

        if self.settings.ensemble == "NVE":

            #Compute all relevant quantities and write them to a file 
            pass
            
        if self.settings.ensemble == "NVE":

            #Compute all relevant quantities and write them to a file 
            pass
            
        if self.settings.ensemble == "NVE":

            #Compute all relevant quantities and write them to a file 
            pass
        

    def writeQuantities(self,labels,quantities):
        """
        Write labels and quantities to txt file. Currently doesnt support timeseries
        Ex:

        label1  label2 ....
        q1      q2    ......
        
        """
        col_width = 12
        with open(f"{self.settings.output_file}.txt", "w") as f:
            f.write(f"Ensemble: {self.settings.ensemble}\n")
            #TODO Add more data to the header
            
            f.write("".join(f"{label:<{col_width}}" for label in labels) + "\n")
            
            f.write("".join(f"{value:<{col_width}.3f}" for value in quantities) + "\n")


        




    #### CALCULATE ALL QUANTITIES

    def computeMSD(self, time=-10, reference=0):
        # TODO Shouldn't this be some mean value of many of the late snapshots in trajectory
        

        r_0 = self.traj[reference].get_positions()  # Å
        r_n = self.traj[time].get_positions()  # Å
        
        msd = np.mean((r_0 - r_n) ** 2)
        logger.debug(f"MSD: {msd} Å^2")
        return msd
    
    def computeSelfDiffusionCoefficient(self):  # Needs constant temperature, for current implementation
        # Find the actual elapsed time
        timestep_list = []
        for i in range(len(self.traj)):
            timestep = i * self.settings.timestep * self.settings.sample_interval
            timestep_list.append([timestep, i])

        if len(timestep_list) > 100:
            msd0 = self.computeMSD(time=timestep_list[50][1])
            msd_final = self.computeMSD(time=timestep_list[-1][1])
            t_0 = timestep_list[50][0]
            t_end = timestep_list[-1][0]
            D = (msd_final-msd0)/(t_end - t_0)
            logger.debug(f"Self-diffusion coefficent:{D}")
        else:
            logger.error("Too small sample size to calcualte self-diffusion coefficient")
            D = None
        
        
        return D * 10**-5 / 6

