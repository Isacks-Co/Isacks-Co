from ase.io.trajectory import Trajectory
from simulationInput import SimulationInput


class QuantityCalculator:

    def __init__(self,settings: SimulationInput, traj: Trajectory):
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




