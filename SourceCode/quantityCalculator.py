from SourceCode.simulationInput import SimulationSettings
from SourceCode.Utils.unitConversions import auToGPascal,specificHeatAuToSI,selfDiffusionCoeffAuToSI, evToJ
from SourceCode.Utils.plotting import secondOrderNumericalDerivative, numericalDerivative

from scipy.constants import physical_constants
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.calculators.emt import EMT
from ase.units import kB
from ase.eos import EquationOfState
from matplotlib import pyplot as plt
from collections import defaultdict

import numpy as np
import logging
from potentialSetUp import Potential

hbar = physical_constants['Planck constant over 2 pi in eV s'][0] * 1e15

logger = logging.getLogger(__name__)

class QuantityCalculator: 
    """
    Object for handling computation of quantities. 
    """
    

    @staticmethod
    def computeCohesiveEnergy(atoms):
        """
        Calculate the cohesive energy per atom. 
        This is done by computing the difference in energy between the separate atoms and the bulk structure 
        
        Unit: ev/atom
        """
        
        number = atoms.get_global_number_of_atoms()
        e_atoms = 0

        for symbol in atoms.get_chemical_symbols():
            atom = Atoms(symbol, positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=False)

            atom.calc = atoms.calc
            e_atoms += atom.get_potential_energy()
        
        e_bulk = _get(atoms, "E_pot")
        e_coh = (e_bulk-e_atoms) / number
        logger.info(f"Cohesive energy: {e_coh} eV")
        return e_coh

    @staticmethod
    def computeSpecificHeatNVT(atoms_sequence): # TODO Make decision on where to put this
        """ 
        Compute specific heat capacity as a time average of the total energy fluctuations 
        over all frames.
        Unit: ev/(amu*K) (amu = atomic mass unit)
        """
        #energy = np.array([atom_frame.get_total_energy() for atom_frame in self.traj])
        #temperature = np.mean([atom_frame.get_temperature() for atom_frame in self.traj])
        energy = np.array([_get(frame, "E_tot") for frame in atoms_sequence])
        temperature = np.mean([_get(frame, "T") for frame in atoms_sequence])

        e_mean = np.mean(energy)
        e_2_mean = np.mean(energy**2)
        prefactor = 1/(kB*temperature**2)
        total_mass_amu = float(sum(atoms_sequence[0].get_masses()))
        specific_heat =  prefactor * (e_2_mean-e_mean**2)/total_mass_amu # Specific heat in ev/amu*K
        logger.debug(f"Cv: {specificHeatAuToSI(specific_heat)} J/(kg*K)")
        return specific_heat
    
    @staticmethod
    def computeSpecificHeatNVE(atoms_sequence): # TODO Make decision on where to put this
        """
        Compute specific heat capacity as a time average of the kinetic energy fluctuations 
        over all frames. 
        Unit: ev/(amu*K) (amu = atomic mass units)
        """
        e_kin = np.array([_get(fr, "E_kin") for fr in atoms_sequence])
        T = float(np.nanmean([_get(fr, "T") for fr in atoms_sequence]))

        e_kin_mean = np.mean(e_kin)
        e_kin_2_mean = np.mean(e_kin**2)
        total_mass_amu = float(sum(atoms_sequence[0].get_masses()))
        specific_heat =  (3*kB/2)*1/(1-(2/(3*(kB*T)**2)*(e_kin_2_mean - e_kin_mean**2)))/total_mass_amu # Should verify
        logger.debug(f"Cv: {specificHeatAuToSI(specific_heat)} J/(kg*K)")
        return specific_heat
    
    @staticmethod
    def computeLatticeConstant(atoms): # Need to test
        #lattice_frames = [atoms.get_cell().cellpar() for atoms in self.traj]
        
        #lattice_mean = np.mean(lattice_frames,axis = 0)

        lattice_mat = atoms.get_cell().cellpar()
        logger.info(f"Lattice constant: {lattice_mat}")
        return lattice_mat
    
    @staticmethod
    def computeInternalPressure(atoms_sequence): # TODO Currently only valid for NVT look into implementaions for other ensembles
        """
        Compute internal pressure using atomic units internally.
        Instantaneous: P = (1/3V) [ 2 N E_kin + sum_i r_i · f_i ]
        Returns a mean over the instantaneous frames in the trajectory
        Unit: ev/Å^3
        """
        internal_pressures_eVA3 = []
        N = len(atoms_sequence[0])

        for atoms in atoms_sequence:
            e_kin_eV = _get(atoms, "E_kin")
            V_A3 = _get(atoms, "V")
            forces_eVA = _get(atoms, "F")
            positions_A = atoms.get_positions()
            sum_rf = np.sum(forces_eVA * positions_A)
            P_eVA3 = (1.0 / (3.0 * V_A3)) * (2.0 * e_kin_eV + sum_rf)
            internal_pressures_eVA3.append(P_eVA3)

        avg_P = np.mean(internal_pressures_eVA3) if internal_pressures_eVA3 else float('nan')
        logger.debug(f"Average internal pressure: {avg_P} eV/Å^3")
        return avg_P
    @staticmethod
    def computeMSD( atoms_sequence,frame, reference=0):
        """
        Compute MSD for a specific frame relatice the first frame in the trajectory. 
        Should not use data from the first 10 points since to close to reference
        Unit: Å^2
        """
        r_0 = atoms_sequence[reference].get_positions()  # Å
        r_n = atoms_sequence[frame].get_positions()  # Å

        msd = np.mean((r_0 - r_n) ** 2)
        logger.debug(f"MSD: {msd} å²")
        return msd
    
    @staticmethod
    def computeSelfDiffusionCoefficient(atoms_sequence,sample_spacing):  # Needs constant temperature, for current implementation. 
        
        """
        Compute self diffusion coefficient from the slope of MSD over a large timeperiod
        Unit: Å^2/fs
        """

        # Find the actual elapsed time
        timestep_list = []
        for i in range(len(atoms_sequence)):
            timestep = i * sample_spacing
            timestep_list.append([timestep, i])

        if len(timestep_list) > 100: # Since early values of MSD are inaccurate
            msd0 = QuantityCalculator.computeMSD(atoms_sequence=atoms_sequence,frame=timestep_list[50][1])
            msd_final = QuantityCalculator.computeMSD(atoms_sequence=atoms_sequence,frame=timestep_list[-1][1])
            t_0 = timestep_list[50][0]
            t_end = timestep_list[-1][0]
            D = (msd_final-msd0)/(t_end - t_0)
            logger.debug(f"Self-diffusion coefficent: {selfDiffusionCoeffAuToSI(D)} m²/s")
        else:
            logger.error("Too small sample size to calculate self-diffusion coefficient")
            D = None

        return D / 6 # Å^2/fs





    #TODO Still need to look at the stuf below this
    

    @staticmethod
    def nearestNeighborsMean( atoms_sequence, start: int, end: int = None): # TODO Look at this and make work without traj
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
            atoms = atoms_sequence[state]
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
        logger.debug(f"Mean value of nearest neighbor : {NN_mean_distance} å")
        return NN_mean_distance
    
    @staticmethod
    def computeLindemannIndex( atoms_sequence, start:int = -25, end:int = 0):
        """Returns the global Lindemann index for the given interval
        (int) start : index for the start of the interval that should be checked
        (int) end  : index for the end of the interval that should be checked
        """
        lindemann_array = []
        for state in range(start, end):
            lindemann_array.append(np.sqrt(QuantityCalculator.computeMSD(atoms_sequence,frame = state)) / QuantityCalculator.nearestNeighborsMean(atoms_sequence,state))
        lindemann = np.mean(lindemann_array)

        logger.debug(f"Global Lindemann index for the intervals [{start}, {end}] : {lindemann}")
        return lindemann



    @staticmethod
    def computeDebyeTemperature(atoms, G, K):
        V_A3 =atoms.get_volume()
        mass_u = float(sum(atoms.get_masses()))
        rho = (mass_u / V_A3)

        transversal_sound_velocity = np.sqrt(G / rho)
        longitudinal_sound_velocity = np.sqrt((K + 4.0 * G / 3.0) / rho)
        sound_velocity = ((1.0 / 3.0) * (1.0 / (longitudinal_sound_velocity ** 3) + 2.0 / (transversal_sound_velocity ** 3))) ** (-1.0 / 3.0)

        N = len(atoms)
        n = (N / V_A3)

        Theta_D = (hbar / kB) * ((6.0 * np.pi ** 2 * n) ** (1.0 / 3.0)) * sound_velocity / 10.18 # NEED TO DO SQRT(ev/u) to fs/Å
        logger.info(f"Debye temperature: {Theta_D} K")
        return Theta_D

    @staticmethod
    def calculateModuli(C_matrix):
        bulk_modulus = (C_matrix[0, 0] + 2 * C_matrix[0, 1]) / 3
        G_shear = (C_matrix[3, 3] + C_matrix[4, 4] + C_matrix[5, 5] + C_matrix[1, 1] - C_matrix[0, 1]) / 5
        youngs_modulus = 9 * bulk_modulus * G_shear / (3 * bulk_modulus + G_shear)
        return bulk_modulus, G_shear, youngs_modulus





    

    @staticmethod
    def computeBulkModulus( stretch_sequence):
        energies = []
        cells = []
        for frame in stretch_sequence:
            energies.append(_get(frame, "E_pot"))
            cells.append(_get(frame, "V"))

        V = np.array(cells)
        E = np.array(energies)
        order = np.argsort(V)
        E, V = E[order], V[order]

        eos = EquationOfState(V, E, eos='birchmurnaghan')
        v0, e0, B0_eVa3 = eos.fit()
        B0_GPa = auToGPascal(B0_eVa3)
        eos.plot('Ag-eos.png')
        return B0_GPa






    @staticmethod
    def calculateCMatrix(strech_sequence):
        
        betas = [[], [], [], [], [], []]
        for frame in strech_sequence:
            # Add the corresponding values to the right beta (strain direction)
            betas[frame.info["beta"]].append([frame.info["strain"], frame.info["stress"]])
        beta_arrays = [np.array(beta, dtype=object) for beta in betas]

        # Create dictionaries, one for each beta, and store the arrays of matrices with epsilon as key
        beta_dicts = [defaultdict(list) for i in range(6)]
        for i in range(6):
            for epsilon, matrix in beta_arrays[i]:
                beta_dicts[i][epsilon].append(matrix)

        averages = []
        for beta in beta_dicts:
            avg_data = []

            for epsilon, matrices in beta.items():
                # Take elementwise average over all the matrices for each epsilon
                stacked = np.stack(matrices)
                avg_matrix = stacked.mean(axis=0)
                avg_data.append(np.array((epsilon, avg_matrix), dtype=object))

            avg_data = sorted(avg_data, key=lambda x: x[0])
            averages.append(avg_data)
        C = np.zeros((6, 6))
        for i in range(6):
            # Line fit epsilon vs sigma to find each c_ij
            epsilons = np.array([x[0] for x in averages[i]], dtype=float)
            for j in range(6):
                sigmas = np.array([x[1][j] for x in averages[i]], dtype=float)
                C[j, i] = np.polyfit(epsilons, sigmas, 1)[0]
        C *= 160.21766208  # Convert to GPa
        return C
    


    @staticmethod
    def _numericalC(stretch_sequence): # TODO Not sure if this will work with my changes
        """
        Calculates the elastic constants C11, C22, C33, C12, C44

        Input
            ---

        Output
            C_from_U: matrix, where C11, C12, C44 = C[0,0], C[0,1], C[3,3]
        """
        betas = [[], [], [], [], [], []]
        # Prefer reconstructing 1D slices from 2D trajectory to avoid stale 1D data
        
        used_2d = False
        try:
            
            tol = 1e-12
            for frame in stretch_sequence:
                info = getattr(frame, 'info', {})
                b1 = info.get('beta1')
                b2 = info.get('beta2')
                if b1 is None or b2 is None:
                    continue
                if int(b1) == int(b2) and np.isclose(float(info.get('strain1', 0.0)), 0.0, atol=tol):                   #Looks at difference between strains when creating 2D data
                    i = int(b1)
                    eps = float(info.get('strain2', 0.0))
                    energy = float(info.get('total_energy', np.nan))
                    stress = np.array(info.get('stress'), dtype=float)
                    betas[i].append([eps, stress, energy])
            used_2d = any(len(b) > 0 for b in betas)
        except Exception as e2:
            logger.info(f"Failed")
        if not used_2d:
            # Fallback to legacy 1D if available
            try:
                stretch_trajectory = Trajectory(self.settings.output_file + "_stretch_data.traj")
                for frame in stretch_trajectory:
                    energy = frame.info.get('total_energy')
                    betas[int(frame.info['beta'])].append([float(frame.info['strain']), np.array(frame.info['stress'],
                                                                                        dtype=float), float(energy)])
                logger.info("Used legacy 1D stretch trajectory to compute elastic constants.")
            except Exception as e:
                logger.info(f"No usable stretch data found (2D or 1D): {e}")
        beta_arrays = [np.array(beta, dtype=object) for beta in betas]
        beta_dicts = [defaultdict(list) for i in range(6)]
        averages = []

        for i in range(6):
            for eps, matrix, energy in beta_arrays[i]:
                beta_dicts[i][eps].append((matrix, energy))

        for beta in beta_dicts:
            avg_data = []

            for eps, matrices_and_energies in beta.items():
                # Take elementwise average for each epsilon
                matrices = np.stack([me[0] for me in matrices_and_energies])
                energies = np.array([me[1] for me in matrices_and_energies], dtype=float)
                avg_matrix = matrices.mean(axis=0)
                avg_energy = np.nanmean(energies)  # ignore possible NaNs
                avg_data.append(np.array((eps, avg_matrix, avg_energy), dtype=object))

            avg_data = sorted(avg_data, key=lambda x: x[0])
            averages.append(avg_data)

        second_deriv = np.zeros((6,6))
        

        twoD_energies = stretch_sequence[-1].info["2D Energies"]
        strains_axis = stretch_sequence[-1].info["Strains axis"]
        number_of_pairs = stretch_sequence[-1].info["Number of pairs"]         # Should usually be 6

        for i in range(int(np.sqrt(number_of_pairs))):
            energy_1 = np.array([x[2] for x in averages[i]], dtype=float)
            for j in range(int(np.sqrt(number_of_pairs))):
                energy_2 = np.array([x[2] for x in averages[j]], dtype=float)
                if i == j:
                    second_deriv[i, j] = secondOrderNumericalDerivative(strains_axis, [energy_1, energy_2])
                    continue
                else:
                    try:
                        second_deriv[i, j] = secondOrderNumericalDerivative(strains_axis, twoD_energies[i][j])
                    except Exception as e:
                        logger.info(f"2D stretch calc failed for ({i},{j}), WHY THE FRICK???!: {e}")
                        second_deriv[i, j] = 0.0

        C_from_U = second_deriv / stretch_sequence[0].get_volume()
        logger.debug(f"C_from_U = \n {C_from_U * auToGPascal(1)} \n")
        logger.info(f" \n C_11 = {auToGPascal(C_from_U[0,0])} \n C_12 = {auToGPascal(C_from_U[0,1])} \n C_44 = {auToGPascal(C_from_U[3,3])}")

        return C_from_U












    def computeDebyeTemperature(self): # TODO Problematic
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
            out = self.calcNumericElasticModuli()  # SI Pa
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



   



    

    


