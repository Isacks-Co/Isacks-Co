<<<<<<< HEAD
from ase.io.trajectory import Trajectory
from simulationInput import SimulationSettings
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
import logging

logger = logging.getLogger(__name__)

class QuantityCalculator:

    def __init__(self,settings: SimulationSettings, traj : Trajectory):
        self.traj = traj
=======
from .simulationInput import SimulationSettings

import numpy as np
from ase.units import kB
from ase.io.trajectory import Trajectory

class QuantityCalculator:

    def __init__(self,settings: SimulationSettings):
        self.traj = Trajectory(f"{settings.output_file}.traj")
>>>>>>> 4f1c190 (Refactored some stuff)
        self.settings = settings


    def getQuantities(self):
        # Compute all general quantities
        # Cohesive Energy
        # Internal pressure
        # MSD
        # Lindemann criteria
        # Self diffusion coefficient
        # Specific Heat
        # Debye Temperature

        if self.settings.ensemble == "NVE":

            #Compute all relevant quantities and write them to a file 
            pass
            
        if self.settings.ensemble == "NVT":
            

            #Compute all relevant quantities and write them to a file 
            pass
            
        if self.settings.ensemble == "NPT":
             
            # Lattice Constant
            # Bulk modulus
            # Compute all relevant quantities and write them to a file 
            pass
        



    def computeSpecificHeatNVT(self): #Cv per atom in a.u
        
        energy = np.array([atom_frame.get_total_energy() for atom_frame in self.traj])
        temperature = np.mean([atom_frame.get_temperature() for atom_frame in self.traj]) # Should we do this or just read from settings?
        N = len(self.read_traj_file[0])
        e_mean = np.mean(energy)
        e_2_mean = np.mean(energy**2)
        prefactor = 1/(kB*temperature**2)
        specific_heat =  prefactor * (e_2_mean-e_mean**2)/N
        return specific_heat
    
    def computeSpecificHeatNVE(self): #Cv per atom in a.u
        e_kin = np.array([atom_frame.get_kinetic_energy() for atom_frame in self.traj])
        T = np.mean([atom_frame.get_temperature() for atom_frame in self.traj])
        e_kin_mean = np.mean(e_kin)
        e_kin_2_mean = np.mean(e_kin**2)

        specific_heat =  (3*kB/2)*1/(1-(2/(3*(kB*T)**2)*(e_kin_2_mean - e_kin_mean**2)))
        return specific_heat
    
    def computeBulkModulus(self):
        volumes = np.array([atom_frame.get_volume() for atom_frame in self.traj[2000:]])
        V_mean = np.mean(volumes)
        V_2_mean = np.mean(volumes**2)
        T = np.mean([atom_frame.get_temperature() for atom_frame in self.traj[2000:]])
        B = kB*T*V_mean/(V_2_mean-V_mean**2)
        print("BULK: ", B)
        return B

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


        






    def computeMSD(self, frame, reference=0):
        
        r_0 = self.traj[reference].get_positions()  # Å
        r_n = self.traj[frame].get_positions()  # Å
        
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
            logger.error("Too small sample size to calculate self-diffusion coefficient")
            D = None

        
        return D * 10**-5 / 6


    def nearestNeighborsMean(self, start: int, end: int = None):
        """Calculate the mean distance of nearest neighbor in the structure for the last ten states of the simulation
        Loop structure: Last ten states -> Each atom -> neighbors to current atom

        (int) start : Index for the start of the interval that should be checked
        (int) end : Index for the end of the interval that should be checked
        """
        if end is None:
            end = start + 1
        INF = 1e9
        NN_list = []

        for state in range(start, end):
            # Load neighbor list for the current state
            atoms = self.traj[state]
            cutoff = natural_cutoffs(atoms)
            neighbor_list = NeighborList(cutoff, bothways=True)
            neighbor_list.update(atoms)

            for current_atom in range(atoms.get_global_number_of_atoms()):
                # Loop over all atoms in and find their nearest neighbor
                indices, offsets = neighbor_list.get_neighbors(current_atom)
                nearest_distance = INF

                # First object seems to be the atom itself, don't loop over it
                for neighbor_index, offset in zip(indices[1:], offsets[1:]):
                    # Create a vector between current_atom and the neighbors in the list, save the shortest distance
                    NN_vector = atoms.positions[neighbor_index] + offset @ atoms.get_cell() - atoms.positions[
                        current_atom]
                    distance = np.sqrt(NN_vector.dot(NN_vector))

                    if distance < nearest_distance:
                        nearest_distance = distance

                if nearest_distance == INF:
                    error_msg = f"Could not calculate NN distance, didn't find any NN for atom {current_atom}"
                    logger.error(error_msg)
                    return

                NN_list.append(nearest_distance)
        NN_mean_distance = np.mean(NN_list)
        logger.debug(f"Mean value of nearest neighbor : {NN_mean_distance}")
        return NN_mean_distance


    def computeLindemannIndex(self, start:int = -25, end:int = 0):
        """Returns the global Lindemann index for the given interval
        (int) start : index for the start of the interval that should be checked
        (int) end  : index for the end of the interval that should be checked
        """
        lindemann_array = []
        for state in range(start, end):
            lindemann_array.append(np.sqrt(self.computeMSD(time = state)) / self.nearestNeighborsMean(state))
        lindemann = np.mean(lindemann_array)

        logger.debug(f"Global Lindemann index for the intervals [{start}, {end}] : {lindemann}")
        return lindemann

