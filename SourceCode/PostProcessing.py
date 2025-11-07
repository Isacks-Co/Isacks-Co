import logging
from simulationInput import SimulationSettings
from ase.io.trajectory import Trajectory
from unitConversions import selfDiffusionCoeffAuToSI,auToGPascal,specificHeatAuToSI
from quantityCalculator import QuantityCalculator as QC


class PostProcessing():
    """
    Class responsible for handling the interface between the simulations and
    the postprocessing task such as visualization and computataition of quantities.
    Typically takes one or multiple trajectories and some flags for what to compute.
    
    """



    def __init__(self, trajectory_path):
        
        self.traj = Trajectory(trajectory_path)
        # Converts AND returns functions instead of values
        self.quantity_map = {
            "D": lambda: selfDiffusionCoeffAuToSI(QC.computeSelfDiffusionCoefficient()),
            "E_coh": lambda: QC.computeCohesiveEnergy(),
            "L_crit": lambda: QC.computeLindemannIndex(),

            # NVE-only
            "Cv_NVE": lambda: specificHeatAuToSI(QC.computeSpecificHeatNVE()),

            # NVT-only
            "P_i": lambda: auToGPascal(self.computeInternalPressure()),
            "Cv_NVT": lambda: specificHeatAuToSI(QC.computeSpecificHeatNVT()),
        }

     





    def computeQuantities(self,quantities):
        """
        Compute all relevant quantities for the given ensemble 
        and write them to a txt file
        """
        # Compute all general quantities

        #MSD = self.computeMSD() # Å  Should we output average over late frames ?
        bulk_modulus = QC.computeBulkModulus("../Outputs/isotropic_stretch.traj")
        self_diffusion_coeff = selfDiffusionCoeffAuToSI(QC.computeSelfDiffusionCoefficient())# m^2/s
        coh_energy = QC.computeCohesiveEnergy() # ev
        internal_pressure = auToGPascal(QC.computeInternalPressure()) #GPa
        lindemann_crit = QC.computeLindemannIndex() # Unitless
        #debye_temperature = self.computeDebyeTemperature()

        #logger.info(self.elastic_properties)

        labels = ["D[m^2/s]","E_coh[eV]","L_crit[1]"] #, "T_D"
        quantities = [self_diffusion_coeff,coh_energy,lindemann_crit] # TODO Maybe nicer way to handle this ? , debye_temperature
        match self.settings.ensemble:
            case "NVE":
                Cv = specificHeatAuToSI(QC.computeSpecificHeatNVE()) # J/K per atom
                labels.append("Cv[J/kgK]")
                quantities.append(Cv)

            case "NVT":
                print("NVT")
                internal_pressure = auToGPascal(self.computeInternalPressure())  # GPa
                bulk_modulus = self.computeBulkModulus("../Outputs/isotropic_stretch.traj")
                Cv = specificHeatAuToSI(QC.computeSpecificHeatNVT()) # J/K per atom
                labels.extend(["P_i[GPa]", "B[GPa]", "Cv[J/kgK]"])
                quantities.extend([internal_pressure, bulk_modulus, Cv])


                C_matrix = self.calculateCMatrix()
                bulk_modulus, g_shear, youngs_modulus = self.calculateModuli(C_matrix)
                labels.extend(["B", "G"])
                quantities.extend([bulk_modulus, g_shear])

            case "NPT":
                pass

        self.writeQuantities(labels,quantities) # Write to txt file


    def writeQuantities(self,labels,quantities):
        """
        Write labels and quantities to txt file. 
        Ex:

        HEADER 
        
        label1  label2  .....
        q1      q2      .....
    
        """
        col_width = 20 # Should work fine with current number of decimals
        with open(f"{self.settings.output_file}.txt", "w") as f:
            #HEADER
            f.write(f"{self.structure_name}\n")
            f.write(f"Ensemble: {self.settings.ensemble}\n")
            #TODO Add more data to the header
            f.write(f"\n")

            #DATA
            f.write("".join(f"{label:<{col_width}}" for label in labels) + "\n")
            f.write("".join(f"{value:<{col_width}.3f}" for value in quantities) + "\n")
if __name__ ==  "__main__":
