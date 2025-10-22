from simulationInput import SimulationSettings
from unitConversions import AuToGPascal,specificHeatAuToSI,selfDiffusionCoeffAuToSI

from scipy.constants import physical_constants
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.calculators.emt import EMT
from ase.units import kB
from ase.eos import EquationOfState

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
        #self.elastic_properties = self._cubicConstantsFromTrajectory()

    def getQuantities(self):
        """
        Compute all relevant quantities for the given ensemble 
        and write them to a txt file
        """
        # Compute all general quantities

        #MSD = self.computeMSD() # Å  Should we output average over late frames ?
        bulk_modulus = self.computeBulkModulus("../Outputs/isotropic_stretch.traj")
        self_diffusion_coeff = selfDiffusionCoeffAuToSI(self.computeSelfDiffusionCoefficient())# m^2/s
        coh_energy = self.computeCohesiveEnergy() # ev
        internal_pressure = AuToGPascal(self.computeInternalPressure()) #GPa
        lindemann_crit = self.computeLindemannIndex() # Unitless
        #debye_temperature = self.computeDebyeTemperature()

        #logger.info(self.elastic_properties)

        labels = ["D","E_coh","P_i","L_crit", "B"] #, "T_D"
        quantities = [self_diffusion_coeff,coh_energy,internal_pressure,lindemann_crit, bulk_modulus] # TODO Maybe nicer way to handle this ? , debye_temperature
        match self.settings.ensemble:
            case "NVE":
                Cv = specificHeatAuToSI(self.computeSpecificHeatNVE()) # J/K per atom
                labels.append("Cv")
                quantities.append(Cv)

            case "NVT":
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
            #e_bulk = frame.get_potential_energy()
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
        #energy = np.array([atom_frame.get_total_energy() for atom_frame in self.traj])
        #temperature = np.mean([atom_frame.get_temperature() for atom_frame in self.traj])
        energy = np.array([_get(frame, "E_tot") for frame in self.traj])
        temperature = np.mean([_get(frame, "T") for frame in self.traj])

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
        #e_kin = np.array([atom_frame.get_kinetic_energy() for atom_frame in self.traj])
        #T = np.mean([atom_frame.get_temperature() for atom_frame in self.traj])
        e_kin = np.array([_get(fr, "E_kin") for fr in self.traj])
        T = float(np.nanmean([_get(fr, "T") for fr in self.traj]))

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
        """
        for atoms in self.traj:
            e_kin_eV = atoms.get_kinetic_energy()
            V_A3 = atoms.get_volume()
            forces_eVA = atoms.get_forces()
            positions_A = atoms.get_positions()
            sum_rf_Eh = np.sum(forces_eVA * positions_A)
            P_eVA3 = (1.0 / (3.0 * V_A3)) * (2.0 * e_kin_eV + sum_rf_Eh)
            internal_pressures_eVA3.append(P_eVA3)
        """
        for atoms in self.traj:
            e_kin_eV = _get(atoms, "E_kin")
            V_A3 = _get(atoms, "V")
            forces_eVA = _get(atoms, "F")
            positions_A = atoms.get_positions()

            if forces_eVA is None:
                pass
            #print("FORCE", forces_eVA)
            #print("POSITIONS", positions_A)
            sum_rf_Eh = np.sum(forces_eVA * positions_A)
            P_eVA3 = (1.0 / (3.0 * V_A3)) * (2.0 * e_kin_eV + sum_rf_Eh)
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

        """if self.settings.temperature < 1:
            # Low-temperature Debye from heat capacity; use a.u. for kB
            C_v = self.computeSpecificHeat()  # [J kg-1 K-1]
            temperature = float(np.mean([atom_frame.get_temperature() for atom_frame in self.traj]))
            N = self.traj[0].get_global_number_of_atoms()
            # Here we keep the classical constant form but ensure SI at the end
            debye = (234 * N * kB * EV_TO_JOULE * temperature ** 3 / C_v) ** (1 / 3)
            logger.info(f"Debye temperature: {debye} K")
            return float(debye)"""
    
        # Θ_D = (ħ/kB) (6π^2 n)^(1/3) vm
        out = self.elastic_properties  # SI Pa
        G = out['G']
        K = out['K']
        EV_PER_A3_TO_GPA = 160.21766208
        logger.info(f"Bulk modulus = {K * EV_PER_A3_TO_GPA} GPa")
        logger.info(f" modulus = {G * EV_PER_A3_TO_GPA} GPa")

        if K < 0 or G < 0:
            error_log = "Negative pressure during calulation of Debye temperature"
            logger.error(error_log)
            raise ValueError(error_log)

        if not (np.isfinite(G) and np.isfinite(K)):
            logger.info("Elastic constants not reliable from traj; Debye temperature cannot be computed.")
            return float('nan')

        #V_A3 = np.mean([fr.get_volume() for fr in self.traj])
        V_A3 = float(np.nanmean([_get(fr, "V") for fr in self.traj]))
        mass_u = float(sum(self.traj[0].get_masses()))
        rho = (mass_u / V_A3)

        transversal_sound_velocity = np.sqrt(G / rho)
        longitudinal_sound_velocity = np.sqrt((K + 4.0 * G / 3.0) / rho)
        sound_velocity = ((1.0 / 3.0) * (1.0 / (longitudinal_sound_velocity ** 3) + 2.0 / (transversal_sound_velocity ** 3))) ** (-1.0 / 3.0)

        logger.info(f"Sound Velociy : {sound_velocity}")
        N = len(self.traj[0])
        n = (N / V_A3)

        logger.info(f"N = {N} : V_A3 = {V_A3}")

        Theta_D = (hbar / kB) * ((6.0 * np.pi ** 2 * n) ** (1.0 / 3.0)) * sound_velocity / 10.18 # NEED TO DO SQRT(ev/u) to fs/Å
        logger.info(f"Debye temperature: {Theta_D} K")
        return Theta_D

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






    def _cubicConstantsFromTrajectory(self, ref=None, tol_abs=1e-8, tol_rel=1e-3):
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
        C_11, C_22, C_33 = [], [], []
        C_12 = []
        C_44, C_55, C_66 = [], [], []

        for frame in stretch_trajectory:
            sig_eVA3 = frame.info["stress"]
            current_measurement = frame.info["measurement"]
            current_stretch_matrix = frame.info["stretch_matrix"]

            match current_measurement:
                case "isotropic":
                    C_11.append(sig_eVA3[0] / (current_stretch_matrix[0][0] - 1))
                    logger.info((sig_eVA3[0] / (current_stretch_matrix[0][0] - 1))  * 160.21766208)
                    C_22.append(sig_eVA3[1] / (current_stretch_matrix[1][1] - 1))
                    C_33.append(sig_eVA3[2] / (current_stretch_matrix[2][2] - 1))
                    C_12.append(sig_eVA3[0] / (current_stretch_matrix[1][1] - 1))

                case "shears_xy":
                    C_66.append(sig_eVA3[5] / (current_stretch_matrix[0][1]))
                case "shears_xz":
                    C_55.append(sig_eVA3[4] / (current_stretch_matrix[0][2]))
                case "shears_yz":
                    C_44.append(sig_eVA3[3] / (current_stretch_matrix[1][2]))
                case _:
                    logger.warning("Didn't recognize current stretch matrix type")

        logger.info(f"Steps : {len(C_11)} {len(C_44)} {len(C_55)} {len(C_66)} ")

        B_bulk = (np.mean(C_11) + 2*np.mean(C_12)) * 160.21766208 / 3
        G_shear = (np.mean(C_44) + np.mean(C_55) + np.mean(C_66) + np.mean(C_11) - np.mean(C_12)) * 160.21766208 / 5
        E_young = 9 * B_bulk * G_shear / (3 * B_bulk + G_shear)


        logger.info(f"C11 : {np.mean(C_11) * 160.21766208} , C22 : {np.mean(C_22) * 160.21766208} , C33 : {np.mean(C_33) * 160.21766208}, C12 : {np.mean(C_12) * 160.21766208}")
        logger.info(f"C44 : {np.mean(C_44) * 160.21766208}, C55 : {np.mean(C_55) * 160.21766208}, C66 : {np.mean(C_66) * 160.21766208}")
        logger.info(f"Bulk modulus B : {B_bulk}")
        logger.info(f"Shear modulus G : {G_shear}")
        logger.info(f"Young modulus E : {E_young}")
        """
            e_xx, e_yy, e_zz = eps[0, 0], eps[1, 1], eps[2, 2]
            e_xy, e_xz, e_yz = eps[0, 1], eps[0, 2], eps[1, 2]
            g_xy, g_xz, g_yz = 2 * e_xy, 2 * e_xz, 2 * e_yz

            s_xx, s_yy, s_zz = sig_eVA3[0], sig_eVA3[1], sig_eVA3[2]
            s_yz, s_xz, s_xy = sig_eVA3[3], sig_eVA3[4], sig_eVA3[5]

            # D pairs
            x_D.extend([e_xx - e_yy, e_yy - e_zz, e_zz - e_xx])
            y_D.extend([s_xx - s_yy, s_yy - s_zz, s_zz - s_xx])

            # Hydrostatic
            tr_eps = e_xx + e_yy + e_zz
            sig_h = (s_xx + s_yy + s_zz) / 3.0
            x_S.append(tr_eps)
            y_S.append(sig_h)

            # Shears
            x_xy.append(g_xy)
            y_xy.append(s_xy)
            x_xz.append(g_xz)
            y_xz.append(s_xz)
            x_yz.append(g_yz)
            y_yz.append(s_yz)

        def slopeFiltered(x, y):
            x = np.asarray(x, dtype=float)
            y = np.asarray(y, dtype=float)
            if x.size < 3:
                return float('nan')
            scale = np.max(np.abs(x))
            mask = np.abs(x) > max(tol_abs, tol_rel * scale)
            x2 = x[mask]
            y2 = y[mask]
            if x2.size < 3:
                return float('nan')
            m = np.dot(x2, y2) / np.dot(x2, x2)
            return float(m)

        D = slopeFiltered(x_D, y_D)
        S_3 = slopeFiltered(x_S, y_S)
        S_3 = S_3 * 3.0 if np.isfinite(S_3) else float('nan')

        C44_xy = slopeFiltered(x_xy, y_xy)
        C44_xz = slopeFiltered(x_xz, y_xz)
        C44_yz = slopeFiltered(x_yz, y_yz)

        C44s = [v for v in [C44_xy, C44_xz, C44_yz] if np.isfinite(v)]
        C44 = float(np.mean(C44s)) if C44s else float('nan')

        if np.isfinite(S_3) and np.isfinite(D):
            C11 = (S_3 + 2.0 * D) / 3.0
            C12 = (S_3 - D) / 3.0
            K = S_3 / 3.0
        else:
            C11 = C12 = K = float('nan')

        # VRH shear modulus
        if np.isfinite(D) and np.isfinite(C44):
            G_V = (D + 3.0 * C44) / 5.0
            denom = 4.0 * C44 + 3.0 * D
            if denom != 0:
                G_R = 5.0 * D * C44 / denom
                G = 0.5 * (G_V + G_R) if np.isfinite(G_R) else G_V
            else:
                G = G_V
        else:
            G = C44

        return {
            'C11': C11, 'C12': C12, 'C44': C44, 'K': K, 'G': G,
            'C44_xy': float(C44_xy) if np.isfinite(C44_xy) else float('nan'),
            'C44_xz': float(C44_xz) if np.isfinite(C44_xz) else float('nan'),
            'C44_yz': float(C44_yz) if np.isfinite(C44_yz) else float('nan')
        }
        """


def _smallStrainFromCells(ref, cur):
    C_ref = ref.cell.array
    C_cur = cur.cell.array
    A_ref = C_ref.T
    # logger.info(f"A_ref: {A_ref} \n")
    A_cur = C_cur.T
    # logger.info(f"A_cur: {A_cur} \n")
    # logger.info(f"INV : {np.linalg.inv(A_ref)} \n")
    F = A_cur @ np.linalg.inv(A_ref)
    # logger.info(f"F: {F} \n")
    eps = 0.5 * (F + F.T) - np.eye(3)
    return eps

def _stressEVA3FromInfo(cur):
    sig_v = cur.info["stress eV/A3"]
    if sig_v is None:
        return None
    return np.asarray(sig_v, dtype=float)  # eV/A3



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
            #1) arrays['F']
            F = frame.arrays.get("F")
            if F is not None:
                return np.asarray(F, float)
            #2) info['F']
            Fi = frame.info.get("F", frame.info.get("forces"))
            if Fi is not None:
                return np.asarray(Fi, float)
            #3) fallback
            try:
                return np.asarray(frame.get_forces(), float)
            except Exception:
                return None


    except Exception:
        return None
