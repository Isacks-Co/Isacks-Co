# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import numpy as np
import pandas as pd
import sys
import json

from Cython.Build.Cache import join_path

from ASEWrappers import AtomicStructure
from Utils.unitConversions import selfDiffusionCoeffAuToSI, auToGPascal, GPascalToAu, specificHeatAuToSI
from ase.io.trajectory import Trajectory
from QuantityCalculator import QuantityCalculator as QC


class PostProcessing():
    """
    Class responsible for handling the interface between the simulations and
    the postprocessing task such as visualization and computataition of quantities.
    Typically takes one or multiple trajectories and some flags for what to compute.
    
    """

    def __init__(self, equil_struct, dataframe, C_matrix, quantities_to_compute):
        self.equil_struct = equil_struct
        self.dataframe = dataframe
        self.C_matrix = C_matrix
        self.time_averages = self.computeTimeAverages().to_frame().T
        self.quantities_to_compute = quantities_to_compute

    def storeQuantities(self):
        self.time_averages = self.time_averages.drop(columns=["time"])
        self.derived_quants = self.computeDerivedQuantities()

        self.writeQuantities(pd.concat([self.time_averages, self.derived_quants], axis=1))

    @classmethod
    def fromFiles(cls, folder, settings_path):
        df = pd.read_fwf(f"{folder}/sampledata.txt", skiprows=1)
        C_matrix = np.load(f"{folder}/cmatrix.npy")
        with open(join_path(folder, settings_path), 'r') as file:
            data = json.load(file)
        quantities_to_compute = data["Compute_quantities"]
        equil_struct = AtomicStructure(Trajectory(f"{folder}/Equil.traj")[-1], data["Simulations_config"]["Supercells"])

        return cls(equil_struct, df, C_matrix, quantities_to_compute)

    def computeDerivedQuantities(self):
        """Computes all the needed quantities that are present in quantities_to_compute."""
        data = {}
        debye_flag = False
        for quantity in self.quantities_to_compute:
            match quantity:
                case "E_coh":
                    pass
                case "Moduli":
                    B, G, E = QC.calculateModuli(C_matrix=self.C_matrix)
                    B = auToGPascal(B)
                    G = auToGPascal(G)
                    E = auToGPascal(E)
                    data["B"] = B
                    data["G"] = G
                    data["E"] = E

                case "Lat_const":
                    data["Lat_const"] = self.equil_struct.lattice_constant

                case "CVT":
                    Cv = specificHeatAuToSI(
                        QC.computeSpecificHeatNVT(self.dataframe["E_tot"], sum(self.equil_struct.masses),
                                                  self.time_averages["T"]))
                    data["Cv"] = Cv

                case "Debye":
                    # Requires the Moduli, do after the loop
                    debye_flag = True

                case "L_crit":
                    L_crit = QC.computeLindemannIndex(self.time_averages["MSD"], self.time_averages["NN"])
                    data["L_crit"] = L_crit

                case "D":
                    D = selfDiffusionCoeffAuToSI(
                        QC.computeSelfDiffusionCoefficient(self.dataframe["MSD"].tolist(), self.dataframe["time"][0]))
                    data["D"] = D

                case "E_coh":
                    data["E_coh"] = self.equil_struct.cohesive_energy(self.dataframe["E_pot"])

            if debye_flag:
                T_D = QC.computeDebyeTemperature(self.time_averages["V"], sum(self.equil_struct.masses),
                                                 len(self.equil_struct), GPascalToAu(data["G"]), GPascalToAu(data["E"]))
                data["T_D"] = T_D

        return pd.DataFrame(data)

    def computeTimeAverages(self):
        return self.dataframe.mean()

    def writeQuantities(self, data: pd.DataFrame):
        """
        Write labels and quantities to txt file. 
        """
        data = data.T

        # Format each row as "Name: Value" with fixed-point numbers
        formatted = data.apply(lambda row: f"{row.name}: {row[0]:.6f}", axis=1)

        # Write everything to CSV
        with open("Quantities.csv", "w") as f:
            for line in formatted:
                f.write(line + "\n")  # each formatted row


if __name__ == "__main__":
    post = PostProcessing.fromFiles(sys.argv[1], sys.argv[2])
    post.storeQuantities()
