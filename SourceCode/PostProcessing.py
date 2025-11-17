import pandas as pd
import numpy as np
from ase.io.trajectory import Trajectory

import sys

from quantityCalculator import QuantityCalculator as QC
from ASEWrappers import AtomicStructure
from Utils.unitConversions import selfDiffusionCoeffAuToSI,auToGPascal,specificHeatAuToSI

class PostProcessing():
    """
    Class responsible for handling the interface between the simulations and
    the postprocessing task such as visualization and computataition of quantities.
    Typically takes one or multiple trajectories and some flags for what to compute.
    
    """
    
    def __init__(self,equil_struct,dataframe, C_matrix):
        self.equil_struct = equil_struct
        self.dataframe = dataframe
        self.C_matrix = C_matrix
        self.time_averages = self.computeTimeAverages().to_frame().T
       
    
    def storeQuantities(self):
        self.time_averages = self.time_averages.drop(columns = ["time","MSD"])
        self.derived_quants = self.computeDerivedQuantities()
        
        self.writeQuantities(pd.concat([self.time_averages,self.derived_quants],axis=1))
        
    @classmethod
    def fromFiles(cls,folder):
      
       
        df = pd.read_fwf(f"{folder}/sampledata.txt",skiprows = 1)
        equil_struct = AtomicStructure(Trajectory(f"{folder}/Equil.traj")[-1])
        C_matrix = np.load(f"{folder}/cmatrix.npy")
        return cls(equil_struct,df,C_matrix)
     
    def computeDerivedQuantities(self):
        #Requires sampledata
        D = selfDiffusionCoeffAuToSI(QC.computeSelfDiffusionCoefficient(self.dataframe["MSD"].tolist(),self.dataframe["time"][1]))
        Cv = specificHeatAuToSI(QC.computeSpecificHeatNVT(self.dataframe["E_tot"],sum(self.equil_struct.masses),self.time_averages["T"]))
        
        #Requires C_Matrix
        B,G,E = QC.calculateModuli(C_matrix=self.C_matrix)
        
        T_D = QC.computeDebyeTemperature(self.time_averages["V"],sum(self.equil_struct.masses),len(self.equil_struct),G,E)
        B = auToGPascal(B)
        G = auToGPascal(G)
        E = auToGPascal(E)
        
        data = {"D":D, "Cv":Cv, "B":B, "G":G, "E":E, "T_D" : T_D}

        return pd.DataFrame(data) 

    def computeTimeAverages(self):
        
        return self.dataframe.mean()





    def writeQuantities(self,data: pd.DataFrame):
        """
        Write labels and quantities to txt file. 
        """
        data  = data.T
     

        # Format each row as "Name: Value" with fixed-point numbers
        formatted = data.apply(lambda row: f"{row.name}: {row[0]:.6f}", axis=1)

        

        # Write everything to CSV
        with open("Quantities.csv", "w") as f:

            for line in formatted:
                f.write(line + "\n")             # each formatted row



if __name__ == "__main__":
    post = PostProcessing.fromFiles(sys.argv[1])
    post.storeQuantities()
    