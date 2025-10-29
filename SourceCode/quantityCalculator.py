from simulationInput import SimulationSettings
from unitConversions import auToGPascal,specificHeatAuToSI,selfDiffusionCoeffAuToSI, evToJ
from Utils import secondOrderNumericalDerivative, load_energy_grid_2d

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
    In:
    settings: A SimulationSettings object containing information regarding the relevan MD run
                such as ensemble, temperature ....
    traj: Trajectory object with the data from each frame in the MD run. Each frame has data such as temperatures, energies and forces. 
    
    """
    def __init__(self,settings: SimulationSettings, traj : Trajectory):
        self.traj = traj
        self.settings = settings
        self.structure_name = self.traj[0].info["comment"]
        self.numeric_elastic_properties = self._numericalC()

    def getQuantities(self):
        """
        Compute all relevant quantities for the given ensemble 
        and write them to a txt file
        """
        # Compute all general quantities

        #MSD = self.computeMSD() # Å  Should we output average over late frames ?
        #bulk_modulus = self.computeBulkModulus("../Outputs/isotropic_stretch.traj")
        self_diffusion_coeff = selfDiffusionCoeffAuToSI(self.computeSelfDiffusionCoefficient())# m^2/s
        coh_energy = self.computeCohesiveEnergy() # ev
        internal_pressure = auToGPascal(self.computeInternalPressure()) #GPa
        lattice_constant = self.computeLatticeConstant()
        lindemann_crit = self.computeLindemannIndex() # Unitless

        elastic_moduli = self.calcNumericElasticModuli()

        labels = ["D[m^2/s]","E_coh[eV]","L_crit[1]"] #, "T_D"
        quantities = [self_diffusion_coeff,coh_energy,lindemann_crit] # TODO Maybe nicer way to handle this ? , debye_temperature
        match self.settings.ensemble:
            case "NVE":
                Cv = specificHeatAuToSI(self.computeSpecificHeatNVE()) # J/K per atom
                labels.append("Cv[J/kgK]")
                quantities.append(Cv)

            case "NVT":
                print("NVT")
                internal_pressure = auToGPascal(self.computeInternalPressure())  # GPa
                #bulk_modulus = self.computeBulkModulus("../Outputs/isotropic_stretch.traj")
                Cv = specificHeatAuToSI(self.computeSpecificHeatNVT()) # J/K per atom
                labels.extend(["P_i[GPa]", "B[GPa]", "Cv[J/kgK]"])
                #quantities.extend([internal_pressure, bulk_modulus, Cv])


                #C_matrix = self.calculateCMatrix()
                #bulk_modulus, g_shear, youngs_modulus = self.calculateModuli(C_matrix)
                labels.append("B")
                labels.append("G")
                labels.append("E")
                #quantities.extend([bulk_modulus, g_shear, youngs_modulus])

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
        calc = Potential().getPotential(self.settings.potential)
        number = self.traj[0].get_global_number_of_atoms()
        e_atoms = 0
        for j in range(len(self.traj[0])):
            atom = Atoms(self.traj[0].get_chemical_symbols()[j], positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=False)

            atom.calc = calc(atom)
            e_atoms += atom.get_potential_energy()

        e_coh_list = []
        for frame in self.traj:
            e_bulk = _get(frame, "E_pot")
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
        energy = np.array([_get(frame, "E_tot") for frame in self.traj])
        temperature = np.mean([_get(frame, "T") for frame in self.traj])

        e_mean = np.mean(energy)
        e_2_mean = np.mean(energy**2)
        prefactor = 1/(kB*temperature**2)
        total_mass_amu = float(sum(self.traj[0].get_masses()))
        specific_heat =  prefactor * (e_2_mean-e_mean**2)/total_mass_amu # Specific heat in ev/amu*K
        logger.debug(f"Cv: {specificHeatAuToSI(specific_heat)} J/(kg*K)")
        return specific_heat

    def computeSpecificHeatNVE(self):
        """
        Compute specific heat capacity as a time average of the kinetic energy fluctuations 
        over all frames. 
        Unit: ev/(amu*K) (amu = atomic mass units)
        """
        e_kin = np.array([_get(fr, "E_kin") for fr in self.traj])
        T = float(np.nanmean([_get(fr, "T") for fr in self.traj]))

        e_kin_mean = np.mean(e_kin)
        e_kin_2_mean = np.mean(e_kin**2)
        total_mass_amu = float(sum(self.traj[0].get_masses()))
        specific_heat =  (3*kB/2)*1/(1-(2/(3*(kB*T)**2)*(e_kin_2_mean - e_kin_mean**2)))/total_mass_amu # Should verify
        logger.debug(f"Cv: {specificHeatAuToSI(specific_heat)} J/(kg*K)")
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

    def computeMSD(self, frame, reference=0):
        """
        Compute MSD for a specific frame relatice the first frame in the trajectory. 
        Should not use data from the first 10 points since to close to reference
        Unit: Å^2
        """
        r_0 = self.traj[reference].get_positions()  # Å
        r_n = self.traj[frame].get_positions()  # Å

        msd = np.mean((r_0 - r_n) ** 2)
        logger.debug(f"MSD: {msd} å²")
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
            logger.debug(f"Self-diffusion coefficent: {selfDiffusionCoeffAuToSI(D)} m²/s")
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
        logger.debug(f"Mean value of nearest neighbor : {NN_mean_distance} å")
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
            out = self.calcElasticModuli()  # SI Pa
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

    def _numericalC(self):
        """
        Calculates the elastic constants C11, C22, C33, C12, C44

        Input
            ---

        Output
            C_from_U: matrix, where C11, C12, C44 = C[0,0], C[0,1], C[3,3]
        """
        betas = [[], [], [], [], [], []]
        # Prefer reconstructing 1D slices from 2D trajectory to avoid stale 1D data
        path2d = self.settings.output_file + "_stretch2D_data.traj"
        used_2d = False
        try:
            traj2d = Trajectory(path2d)
            tol = 1e-12
            for frame in traj2d:
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
            logger.info(f"Failed to read 2D stretch trajectory at {path2d}: {e2}")
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
        path2d = self.settings.output_file + "_stretch2D_data.traj"
        for i in range(6):
            epsilons = np.array([x[0] for x in averages[i]], dtype=float)
            energy_1 = np.array([x[2] for x in averages[i]], dtype=float)
            for j in range(6):
                energy_2 = np.array([x[2] for x in averages[j]], dtype=float)
                if len(epsilons) > 1:
                    # Prefer 2D mixed-partial when available
                    if i != j:
                        try:
                            eps_grid, U_grid = load_energy_grid_2d(path2d, i, j)
                        except Exception as e:
                            logger.info(f"Failed loading 2D grid for pair ({i},{j}): {e}")
                            eps_grid, U_grid = None, None
                        if eps_grid is not None and not isinstance(eps_grid, tuple):
                            try:
                                second_deriv[i, j] = secondOrderNumericalDerivative(eps_grid, U_grid)
                            except Exception as e:
                                logger.info(f"2D mixed partial calc failed for ({i},{j}), falling back: {e}")
                                # fall back to axis-only approximation if shapes allow; otherwise 0.0
                                if len(energy_2) == len(epsilons) and len(epsilons) > 1:
                                    second_deriv[i, j] = secondOrderNumericalDerivative(epsilons, [energy_1, energy_2])
                                else:
                                    second_deriv[i, j] = 0.0
                        else:
                            # No 2D grid available; only safe value for mixed partial is 0.0 unless we have two traces of equal length
                            if len(energy_2) == len(epsilons) and len(epsilons) > 1:
                                second_deriv[i, j] = secondOrderNumericalDerivative(epsilons, [energy_1, energy_2])
                            else:
                                second_deriv[i, j] = 0.0
                    else:
                        # Diagonals: use 1D second derivative path (handled in Utils)
                        second_deriv[i, j] = secondOrderNumericalDerivative(epsilons, [energy_1, energy_2])

        C_from_U = second_deriv / self.traj[0].get_volume()
        logger.debug(f"C_from_U = \n {C_from_U * auToGPascal(1)} \n")
        logger.debug(f" \n C_11 = {auToGPascal(C_from_U[0,0])} \n C_12 = {auToGPascal(C_from_U[0,1])} \n C_44 = {auToGPascal(C_from_U[3,3])}")

        return C_from_U


    def calcNumericElasticModuli(self):
        C = self.numeric_elastic_properties
        K = (C[0,0] + 2 * C[0,1]) / 3
        G = (3 * C[3,3] + C[0,0] - C[0,1]) / 5
        E = 9 * K * G / (3 * K + G)
        logger.debug(f" \n K = {auToGPascal(K)} GPa \n G = {auToGPascal(G)} GPa \n E = {auToGPascal(E)} GPa")
        return {"K": K, "G": G, "E": E}


    def calculateCMatrix(self):
        stretch_trajectory = Trajectory(self.settings.output_file + "_stretch_data.traj")
        betas = [[], [], [], [], [], []]
        for frame in stretch_trajectory:
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

    def calculateModuli(self, C_matrix):
        bulk_modulus = (C_matrix[0, 0] + 2 * C_matrix[0, 1]) / 3
        G_shear = (C_matrix[3, 3] + C_matrix[4, 4] + C_matrix[5, 5] + C_matrix[1, 1] - C_matrix[0, 1]) / 5
        youngs_modulus = 9 * bulk_modulus * G_shear / (3 * bulk_modulus + G_shear)
        return bulk_modulus, G_shear, youngs_modulus

    def computeBulkModulus(self, stretched_traj):
        import matplotlib.pyplot as plt
        stretch_trajectory = Trajectory(stretched_traj)
        energies = []
        cells = []
        for frame in stretch_trajectory:
            energies.append(_get(frame, "E_pot"))
            cells.append(_get(frame, "V"))

        V = np.array(cells)
        E = np.array(energies)
        order = np.argsort(V)
        E, V = E[order], V[order]

        """
        plt.figure(figsize=(6, 4))
        plt.scatter(V, E, s=28)
        plt.xlabel('volume [Å³]')
        plt.ylabel('energy [eV]')
        plt.title('E(V) – stretch trajectory')
        plt.tight_layout()
        plt.savefig('e_vs_v_points.png', dpi=200)
        plt.close()
        """
        eos = EquationOfState(V, E, eos='birchmurnaghan')
        v0, e0, B0_eVa3 = eos.fit()
        B0_GPa = B0_eVa3 * 160.21766208
        eos.plot('Ag-eos.png')
        return B0_GPa


def _get(frame, name):
    """Read from info/arrays first (traj-safe), else fall back to get_* (requires calc)."""
    try:
        if name == "E_pot":
            E_pot = frame.info.get("E_pot")
            return float(E_pot) if E_pot is not None else frame.get_potential_energy()

        if name == "E_kin":
            E_kin = frame.info.get("E_kin")
            return float(E_kin) if E_kin is not None else frame.get_kinetic_energy()

        if name == "E_tot":
            E_tot = frame.info.get("E_tot")
            return float(E_tot) if E_tot is not None else frame.get_total_energy()

        if name == "V":
            v = frame.info.get("V")
            return float(v) if v is not None else frame.get_volume()

        if name == "T":
            T = frame.info.get("T")
            return float(T) if T is not None else frame.get_temperature()

        if name == "F":
            # 1) arrays['F']
            F = frame.arrays.get("F")
            if F is not None:
                return np.asarray(F, float)
            # 2) info['F']
            Fi = frame.info.get("F", frame.info.get("forces"))
            if Fi is not None:
                return np.asarray(Fi, float)
            # 3) fallback
            try:
                return np.asarray(frame.get_forces(), float)
            except Exception:
                return None


    except Exception:
        return None