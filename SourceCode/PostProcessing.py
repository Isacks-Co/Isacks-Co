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
    Post-process simulation outputs to compute averaged and derived quantities.

    The typical workflow is:
    1. Load sample time series (e.g., ``sampledata.txt``) and optional matrices
       (e.g., ``cmatrix.npy``).
    2. Compute time averages over sampled quantities.
    3. Compute derived quantities requested in the settings.
    4. Write results to ``Quantities.csv``.

    Parameters
    ----------
    equil_struct : AtomicStructure
        Equilibrated structure used for structural quantities (e.g., lattice constant).
    dataframe : pandas.DataFrame
        Sampled time series containing columns such as energy, MSD, volume, etc.
    C_matrix : numpy.ndarray
        Elastic stiffness matrix (typically 6×6 in Voigt notation).
    quantities_to_compute : list[str]
        List of quantity identifiers that control what is computed.
    """

    def __init__(self, equil_struct, dataframe, C_matrix, quantities_to_compute):
        self.equil_struct = equil_struct
        self.dataframe = dataframe
        self.C_matrix = C_matrix
        self.time_averages = self.computeTimeAverages().to_frame().T
        # Compute time averages immediately; stored as a 1-row DataFrame for easy concatenation.
        self.quantities_to_compute = quantities_to_compute

    def storeQuantities(self):
        """
        Compute and write requested quantities.

        Drops the ``time`` column from the time-average table, computes derived
        quantities based on `quantities_to_compute`, concatenates results, and
        writes them to ``Quantities.csv``.
        """
        self.time_averages = self.time_averages.drop(columns=["time"])
        self.derived_quants = self.computeDerivedQuantities()

        self.writeQuantities(pd.concat([self.time_averages, self.derived_quants], axis=1))

    @classmethod
    def fromFiles(cls, folder, settings_path):
        """
        Construct a PostProcessing instance from files on disk.

        Parameters
        ----------
        folder : str
            Folder containing MD output files (e.g., ``sampledata.txt``, ``cmatrix.npy``,
            and trajectory files).
        settings_path : str
            Path to the JSON settings file relative to `folder`.

        Returns
        -------
        PostProcessing
            Initialized instance ready for computation and output.
        """

        df = pd.read_fwf(f"{folder}/sampledata.txt", skiprows=1)
        C_matrix = np.load(f"{folder}/cmatrix.npy")
        with open(join_path(folder, settings_path), 'r') as file:
            data = json.load(file)
        quantities_to_compute = data["Compute_quantities"]
        equil_struct = AtomicStructure(Trajectory(f"{folder}/Equil.traj")[-1], data["Simulations_config"]["Supercells"])
        equil_struct.potential = getPotential(data)
        return cls(equil_struct, df, C_matrix, quantities_to_compute)

    def computeDerivedQuantities(self):
        """
        Compute derived quantities requested by `quantities_to_compute`.

        Returns
        -------
        pandas.DataFrame
            Single-row dataframe containing derived quantities as columns.

        Notes
        -----
        Some quantities depend on others (e.g., Debye temperature depends on
        elastic moduli). These dependencies are handled via internal flags.
        """
        data = {}
        debye_flag = False
        for quantity in self.quantities_to_compute:
            match quantity:

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
                    data["E_coh"] = self.equil_struct.cohesive_energy()

            if debye_flag:
                T_D = QC.computeDebyeTemperature(self.time_averages["V"], sum(self.equil_struct.masses),
                                                 len(self.equil_struct), GPascalToAu(data["G"]), GPascalToAu(data["E"]))
                data["T_D"] = T_D

        return pd.DataFrame(data)

    def computeTimeAverages(self):
        """
        Compute time-averaged values for all columns in the sample dataframe.

        Returns
        -------
        pandas.Series
            Mean of each column.
        """
        return self.dataframe.mean()

    def writeQuantities(self, data: pd.DataFrame):
        """
        Write computed quantities to ``Quantities.csv``.

        Parameters
        ----------
        data : pandas.DataFrame
            DataFrame containing one row of values (columns are quantity names).

        Returns
        -------
        None

        Notes
        -----
        The output format is one line per quantity:

            ``<Name>: <Value>``

        Values are written using fixed-point formatting with 6 decimals.
        """
        data = data.T

        # Format each row as "Name: Value" with fixed-point numbers
        formatted = data.apply(lambda row: f"{row.name}: {row[0]:.6f}", axis=1)

        # Write everything to CSV
        with open("Quantities.csv", "w") as f:
            for line in formatted:
                f.write(line + "\n")  # each formatted row

def getPotential(settings):
    from ASEWrappers import LennardJonesPotential, EMTPotential, MACEPotential
    match settings["Simulations_config"]["Potential"]["Kind"]:

        case "LJ":
            from ase.io import read
            import glob
            from Utils import LJParams

            atoms = read(glob.glob("../SetupFiles/atomic_structure*")[0])

            atomic_num = [(atoms.get_atomic_numbers()[0])]
            atomic_symbols = atoms.get_chemical_symbols()
            epsilon = None
            sigma = None
            rc = None
            ro = None
            if "Material" in settings["Simulations_config"]["Potential"]["Parameters"]:
                material = settings["Simulations_config"]["Potential"]["Parameters"]["Material"]
            else:
                material = sorted(set(atomic_symbols))[0]

            if "epsilon_eV" in settings["Simulations_config"]["Potential"]["Parameters"]:
                epsilon = settings["Simulations_config"]["Potential"]["Parameters"]["epsilon_eV"]
            if "sigma" in settings["Simulations_config"]["Potential"]["Parameters"]:
                sigma = settings["Simulations_config"]["Potential"]["Parameters"]["sigma"]
            if "RC" in settings["Simulations_config"]["Potential"]["Parameters"]:
                rc = settings["Simulations_config"]["Potential"]["Parameters"]["RC"]
            if "RO" in settings["Simulations_config"]["Potential"]["Parameters"]:
                ro = settings["Simulations_config"]["Potential"]["Parameters"]["RO"]

            lj_params = LJParams(material=material, epsilon_eV=epsilon, sigma_A=sigma, rc_A=rc, ro_A=ro)
            return LennardJonesPotential(atomic_numbers=atomic_num, epsilons=[lj_params["epsilon_eV"]],
                                         sigmas=[lj_params["sigma_A"]], rc=lj_params["rc_A"])
        case "EMT":
            return EMTPotential()

        case "MACE":
            from ASEWrappers import MACEPotential
            return MACEPotential(model_path=settings["Simulations_config"]["Potential"]["Parameters"]["Path"])

if __name__ == "__main__":
    post = PostProcessing.fromFiles(sys.argv[1], sys.argv[2])
    post.storeQuantities()
