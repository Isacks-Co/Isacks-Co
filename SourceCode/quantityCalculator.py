from simulationInput import SimulationSettings
from unitConversions import auToGPascal,specificHeatAuToSI,selfDiffusionCoeffAuToSI, a1ToM1, a2ToM2, a3ToM3, atomicMassTokg, auToPascal, evToJ
from Utils import plot, fitter

from scipy.constants import physical_constants
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.calculators.emt import EMT
from ase.units import kB
from matplotlib import pyplot as plt

import numpy as np
import logging

hbar = physical_constants['Planck constant over 2 pi in eV s'][0] * 1e15

logger = logging.getLogger(__name__)

class QuantityCalculator:
    """
    Object for handling computation of quantities. 
    In:
    settings: A SimulationSettings object containing information regarding the relevan MD run
                such as ensemble, temperature ....
    traj: Trajectory object with the data from each frame in the MD run. Each frame has data such as temperatures, energies and forces. 
    
    """
    def __init__(self,settings: SimulationSettings, traj : Trajectory):
        self.traj = traj
        self.settings = settings
        self.structure_name = self.traj[0].info["comment"]
        self.elastic_properties = self._C()
        logger.info(self.elastic_properties)

    def getQuantities(self):
        """
        Compute all relevant quantities for the given ensemble 
        and write them to a txt file
        """
        # Compute all general quantities

        #MSD = self.computeMSD() # Å  Should we output average over late frames ?
        self_diffusion_coeff = selfDiffusionCoeffAuToSI(self.computeSelfDiffusionCoefficient())# m^2/s
        coh_energy = self.computeCohesiveEnergy() # ev
        lattice_constant = self.computeLatticeConstant()
        internal_pressure = auToGPascal(self.computeInternalPressure()) #GPa
        lindemann_crit = self.computeLindemannIndex() # Unitless
        debye_temperature = self.computeDebyeTemperature()

        labels = ["D","E_coh","P_i","L_crit", "T_D"]
        quantities = [self_diffusion_coeff,coh_energy,internal_pressure,lindemann_crit, debye_temperature] # TODO Maybe nicer way to handle this ?
        match self.settings.ensemble:
            case "NVE":
                Cv = specificHeatAuToSI(self.computeSpecificHeatNVE()) # J/K per atom
                labels.append("Cv")
                quantities.append(Cv)

            case "NVT":
                print("NVT")
                Cv = specificHeatAuToSI(self.computeSpecificHeatNVT()) # J/K per atom
                labels.append("Cv")
                quantities.append(Cv)

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

    def computeCohesiveEnergy(self):
        """
        Calculate the cohesive energy per atom. 
        This is done by computing the difference in energy between the separate atoms and the bulk structure 
        The return value is the mean of this over all snapshots
        Unit: ev/atom
        """

        calc = EMT() #TODO NEED TO BE SAME AS ACTUAL MD RUN
        number = self.traj[0].get_global_number_of_atoms()
        e_atoms = 0
        for j in range(len(self.traj[0])):
            atom = Atoms(self.traj[0].get_chemical_symbols()[j], positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=False)
            atom.calc = calc
            e_atoms += atom.get_potential_energy()
        e_coh_list = []
        for frame in self.traj:
            e_bulk = frame.get_potential_energy()
            e_coh_list.append((e_atoms - e_bulk)/number)

        e_coh_mean = np.mean(e_coh_list)
        logger.debug(f"Cohesive energy: {e_coh_mean} eV")
        return e_coh_mean

    def computeSpecificHeatNVT(self):
        """
        Compute specific heat capacity as a time average of the total energy fluctuations 
        over all frames.
        Unit: ev/(amu*K) (amu = atomic mass unit)
        """
        energy = np.array([atom_frame.get_total_energy() for atom_frame in self.traj])
        temperature = np.mean([atom_frame.get_temperature() for atom_frame in self.traj])

        e_mean = np.mean(energy)
        e_2_mean = np.mean(energy**2)
        prefactor = 1/(kB*temperature**2)
        total_mass_amu = float(sum(self.traj[0].get_masses()))
        specific_heat =  prefactor * (e_2_mean-e_mean**2)/total_mass_amu # Specific heat in ev/amu*K
        logger.debug(f"Cv: {specific_heat} eV/(amu*K)")
        return specific_heat

    def computeSpecificHeatNVE(self):
        """
        Compute specific heat capacity as a time average of the kinetic energy fluctuations 
        over all frames. 
        Unit: ev/(amu*K) (amu = atomic mass units)
        """
        e_kin = np.array([atom_frame.get_kinetic_energy() for atom_frame in self.traj])
        T = np.mean([atom_frame.get_temperature() for atom_frame in self.traj])
        e_kin_mean = np.mean(e_kin)
        e_kin_2_mean = np.mean(e_kin**2)
        total_mass_amu = float(sum(self.traj[0].get_masses()))
        specific_heat =  (3*kB/2)*1/(1-(2/(3*(kB*T)**2)*(e_kin_2_mean - e_kin_mean**2)))/total_mass_amu # Should verify
        logger.debug(f"Cv: {specific_heat} eV/(amu*K)")
        return specific_heat


    def computeLatticeConstant(self):# Need to test
        lattice_frames = [atoms.get_cell().cellpar() for atoms in self.traj]
        lattice_mean = np.mean(lattice_frames,axis = 0)
        lattice_mean[:3] /= np.array(self.settings.supercells)
        logger.debug(f"Lattice constant: {lattice_mean}")
        return lattice_mean


    def computeInternalPressure(self):
        """
        Compute internal pressure using atomic units internally.
        Instantaneous: P = (1/3V) [ 2 N E_kin + sum_i r_i · f_i ]
        Returns a mean over the instantaneous frames in the trajectory
        Unit: ev/Å^3
        """
        internal_pressures_eVA3 = []
        N = len(self.traj[0])
        for atoms in self.traj:
            e_kin_eV = atoms.get_kinetic_energy()
            V_A3 = atoms.get_volume()
            forces_eVA = atoms.get_forces()
            positions_A = atoms.get_positions()
            sum_rf = np.sum(forces_eVA * positions_A)
            P_eVA3 = (1.0 / (3.0 * V_A3)) * (2.0 * e_kin_eV + sum_rf)
            internal_pressures_eVA3.append(P_eVA3)

        avg_P = np.mean(internal_pressures_eVA3) if internal_pressures_eVA3 else float('nan')
        logger.debug(f"Average internal pressure: {avg_P} eV/Å^3; {avg_P * 160.21766208} GPa; {avg_P * 160.21766208e9} Pa")
        return avg_P


    def computeMSD(self, frame, reference=0):
        """
        Compute MSD for a specific frame relatice the first frame in the trajectory. 
        Should not use data from the first 10 points since to close to reference
        Unit: Å^2
        """
        r_0 = self.traj[reference].get_positions()  # Å
        r_n = self.traj[frame].get_positions()  # Å

        msd = np.mean((r_0 - r_n) ** 2)
        logger.debug(f"MSD: {msd} Å^2")
        return msd


    def computeSelfDiffusionCoefficient(self):  # Needs constant temperature, for current implementation.
        """
        Compute self diffusion coefficient from the slope of MSD over a large timeperiod
        Unit: Å^2/fs
        """

        # Find the actual elapsed time
        timestep_list = []
        for i in range(len(self.traj)):
            timestep = i * self.settings.timestep * self.settings.sample_interval
            timestep_list.append([timestep, i])

        if len(timestep_list) > 100: # Since early values of MSD are inaccurate
            msd0 = self.computeMSD(frame=timestep_list[50][1])
            msd_final = self.computeMSD(frame=timestep_list[-1][1])
            t_0 = timestep_list[50][0]
            t_end = timestep_list[-1][0]
            D = (msd_final-msd0)/(t_end - t_0)
            logger.debug(f"Self-diffusion coefficent:{D}")
        else:
            logger.error("Too small sample size to calculate self-diffusion coefficient")
            D = None

        return D / 6 # Å^2/fs


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
            lindemann_array.append(np.sqrt(self.computeMSD(frame = state)) / self.nearestNeighborsMean(state))
        lindemann = np.mean(lindemann_array)

        logger.debug(f"Global Lindemann index for the intervals [{start}, {end}] : {lindemann}")
        return lindemann


    def computeDebyeTemperature(self):
        """
        Compute Debye temperature using atomic units internally.
        Returns Theta_D in Kelvin (SI).
        """

        if self.settings.temperature < 50:
            # Low-temperature Debye from heat capacity; use a.u. for kB
            if self.settings.ensemble == "NVT":
                C_v = self.computeSpecificHeatNVT()  # [J kg-1 K-1]
            elif self.settings.ensemble == "NVE":
                C_v = self.computeSpecificHeatNVE()  # [J kg-1 K-1]

            temperature = float(np.mean([atom_frame.get_temperature() for atom_frame in self.traj]))
            N = self.traj[0].get_global_number_of_atoms()
            # Here we keep the classical constant form but ensure SI at the end
            debye = (234 * N * kB * evToJ(1) * temperature ** 3 / C_v) ** (1 / 3)
            logger.debug(f"Debye temperature: {debye} K")
            return float(debye)

        else:
            # Θ_D = (ħ/kB) (6π^2 n)^(1/3) vm
            out = self.elastic_properties  # SI Pa
            G = out['G']
            K = out['K']

            if K < 0 or G < 0:
                error_log = "Negative pressure during calulation of Debye temperature"
                logger.error(error_log)
                raise ValueError(error_log)

            if not (np.isfinite(G) and np.isfinite(K)):
                logger.info("Elastic constants not reliable from traj; Debye temperature cannot be computed.")
                return float('nan')

            V_A3 = np.mean([fr.get_volume() for fr in self.traj])
            mass_u = float(sum(self.traj[0].get_masses()))
            rho = (mass_u / V_A3)

            transversal_sound_velocity = np.sqrt(G / rho)
            longitudinal_sound_velocity = np.sqrt((K + 4.0 * G / 3.0) / rho)
            sound_velocity = ((1.0 / 3.0) * (1.0 / (longitudinal_sound_velocity ** 3) + 2.0 / (transversal_sound_velocity ** 3))) ** (-1.0 / 3.0)

            N = len(self.traj[0])
            n = (N / V_A3)

        Theta_D = (hbar / kB) * ((6.0 * np.pi ** 2 * n) ** (1.0 / 3.0)) * sound_velocity / 10.18 # NEED TO DO SQRT(ev/u) to fs/Å
        logger.debug(f"Debye temperature: {Theta_D} K")
        return Theta_D

    def _C(self):
        stretch_trajectory = Trajectory(self.settings.output_file + "_stretch_data.traj")
        xx_dir = []
        xy_dir = []
        for frame in stretch_trajectory:
            if frame.info["beta"] == 0:
                xx_dir.append([frame.info["strain"], frame.info["stress"][0]])
                xy_dir.append([frame.info["strain"], frame.info["stress"][1]])

        xx_dir = np.asarray(xx_dir)
        xy_dir = np.asarray(xy_dir)

        unique_keys, inverse = np.unique(xx_dir[:, 0], return_inverse=True)

        # Compute the mean of the second column for each group
        means = np.bincount(inverse, weights=xx_dir[:, 1]) / np.bincount(inverse)

        # Combine into a result array
        result = np.column_stack((unique_keys, means))
        xx_dir = result
        unique_keys, inverse = np.unique(xy_dir[:, 0], return_inverse=True)

        # Compute the mean of the second column for each group
        means = np.bincount(inverse, weights=xy_dir[:, 1]) / np.bincount(inverse)

        # Combine into a result array
        result = np.column_stack((unique_keys, means))
        xy_dir = result

        C_11 = np.polyfit(xx_dir[:, 0], xx_dir[:,1], 1)[0] * 160.21766208
        C_12 = np.polyfit(xy_dir[:, 0], xy_dir[:,1], 1)[0] * 160.21766208
        plt.scatter(xx_dir[:,0], xx_dir[:,1])
        plt.show()

        logger.info(f"Bulk : {(C_11 + 2*C_12)/3}")
        return C_11, C_12





    def _cubicConstantsFromTrajectory(self, ref=None):
        """
        Estimate C11, C12, C44, K, G using only cell geometry and info['stress'] saved in traj.
        Internally uses atomic units (stress in Eh/a0^3), and returns SI Pascals.
        Voigt notation: 1 -> xx
                        2 -> yy
                        3 -> zz
                        4 -> yz
                        5 -> xz
                        6 -> xy
        """
        stretch_trajectory = Trajectory(self.settings.output_file + "_stretch_data.traj")

        # Collect arrays
        C_11_sigma, C_11_epsilon, C_22, C_33 = [], [], [], []
        C_12_sigma, C_12_epsilon = [], []
        C_44, C_55, C_66 = [], [], []

        for frame in stretch_trajectory:
            sig_eVA3 = frame.info["stress"]
            current_measurement = frame.info["measurement"]
            current_stretch_matrix = frame.info["stretch_matrix"]
            current_energy = frame.info["potential_energy"]

            match current_measurement:
                case "stretch_xx":
                    C_11_sigma.append(sig_eVA3[0])
                    C_11_epsilon.append(current_stretch_matrix[0][0] - 1)
                    C_12_sigma.append(sig_eVA3[1])
                    C_12_epsilon.append(current_stretch_matrix[0][0] - 1)
                    # C_11.append(sig_eVA3[0] / (current_stretch_matrix[0][0] - 1))
                    # C_12.append(sig_eVA3[1] / (current_stretch_matrix[0][0] - 1))
                    # logger.info((sig_eVA3[0] / (current_stretch_matrix[0][0] - 1)) * 160.21766208)
                    # logger.info(f"Energy : {current_energy} , Sigma : {sig_eVA3[0]} , Epsilon : {current_stretch_matrix[0][0] - 1}")

                case "shears_xy":
                    C_66.append(sig_eVA3[5] / (current_stretch_matrix[0][1]))
                case "shears_xz":
                    C_55.append(sig_eVA3[4] / (current_stretch_matrix[0][2]))
                case "shears_yz":
                    C_44.append(sig_eVA3[3] / (current_stretch_matrix[1][2]))
                case _:
                    logger.warning("Didn't recognize current stretch matrix type")

        C_11 = np.polyfit(C_11_epsilon, C_11_sigma, 1)[0]
        C_12 = np.polyfit( C_12_epsilon,C_12_sigma, 1)[0]

        plt.scatter(C_11_epsilon, C_11_sigma, c='r', marker='o', label='C11')
        plt.scatter(C_12_epsilon, C_12_sigma, c='b', marker='o', label='C12')
        plt.plot(C_11_epsilon, C_11 * np.asarray(C_11_epsilon), c='g', label='C11_fit')
        plt.plot(C_12_epsilon, C_12 * np.asarray(C_12_epsilon), c='y', label='C12_fit')
        plt.legend(loc='best')
        plt.show()

        B_bulk = (C_11 + 2*C_12) * 160.21766208 / 3
        G_shear = (np.mean(C_44) + np.mean(C_55) + np.mean(C_66) + np.mean(C_11) - np.mean(C_12)) * 160.21766208 / 5
        E_young = 9 * B_bulk * G_shear / (3 * B_bulk + G_shear)


        logger.info(f"C11 : {C_11 * 160.21766208} , C12 : {C_12 * 160.21766208}")
        logger.info(f"C44 : {np.mean(C_44) * 160.21766208}, C55 : {np.mean(C_55) * 160.21766208}, C66 : {np.mean(C_66) * 160.21766208}")
        logger.info(f"Bulk modulus B : {B_bulk}")
        logger.info(f"Shear modulus G : {G_shear}")
        logger.info(f"Young modulus E : {E_young}")


def _smallStrainFromCells(ref, cur):
    C_ref = ref.cell.array
    C_cur = cur.cell.array
    A_ref = C_ref.T
    A_cur = C_cur.T
    F = A_cur @ np.linalg.inv(A_ref)
    eps = 0.5 * (F + F.T) - np.eye(3)
    return eps


def numericalDerivative(x, y):
    """
    Compute the numerical derivative dy/dx given lists of x and y values.
    Uses central difference for interior points and forward/backward for endpoints.

    Parameters:
        x (list or array): x-values (must be increasing and same length as y)
        y (list or array): y-values

    Returns:
        list: derivative values (same length as x)
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    if len(x) < 2:
        raise ValueError("At least two points are required to compute a derivative.")

    dydx = [0.0] * len(x)

    # Forward difference for first point
    dydx[0] = (y[1] - y[0]) / (x[1] - x[0])

    # Central difference for interior points
    for i in range(1, len(x) - 1):
        if x[i + 1] != x[i - 1]:
            dx = x[i + 1] - x[i - 1]
            dy = y[i + 1] - y[i - 1]
            #logger.debug(f"\n [{i}] ----- dx = {dx} \n [{i}] ----- dy = {dy} \n")
            dydx[i] = dy / dx
        else:
            dydx[i] = (y[i + 1] - y[i - 1]) / (2 * (x[i + 1] - x[i]))
            #logger.info(f"Warning: x values are not unique at index {i}, using central difference.")

    # Backward difference for last point
    dydx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

    return dydx


def sortRealignAndFilter(x_list, y_list, filter=True, resolution=None):
    x_sorted = sorted(x_list)
    y_realigned = []
    for x in x_sorted:
        y_realigned.append(y_list[x_list.index(x)])

    poly, error, x_fit, y_fit = fitter(x_sorted, y_realigned, deg=20, resolution=resolution)
    error_mean = np.mean(list(map(abs, error)))

    if filter:
        to_remove = []
        for j, (quantity, index) in enumerate(zip(y_realigned, x_sorted)):
            if (error[j] / error_mean) >= 2.0:
                logger.debug(f"error ratio at step {j} : {error[j] / error_mean})")
                to_remove.append(j)

        for i in sorted(to_remove, reverse=True):
            print(f"Removing index {i} -> {y_realigned[i]}")
            del y_realigned[i]
            del x_sorted[i]
    return x_sorted, y_realigned