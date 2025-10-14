import json
import re
import math
from pathlib import Path

import numpy as np
from SourceCode.logger import logger_setup
from ase import Atoms
from ase.build import bulk
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.md.analysis import DiffusionCoefficient
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.units import Bohr, kB, Hartree
from ase.visualize import view
from scipy.constants import pi, k, physical_constants
from ase.calculators.emt import EMT


# Common conversions
EV_PER_A3_TO_GPA = 160.21766208
EV_PER_A3_TO_PA = 160.21766208e9  # Pa per (eV/Å^3)
AMU_TO_KG = 1.66053906660e-27
EV_TO_JOULE = 1.602176634e-19
JOULE_TO_EV = 1.0 / EV_TO_JOULE
AVOGRADO = 6.02214076e23
A_TO_M = 1e-10
hbar = physical_constants['Planck constant over 2 pi in eV s'][0]

logger = logger_setup()

class PostProcessing:
    """
        Uses the log files generates from the Molecular Dynamics class to vizualize and analyze the data
        input trajectory file as trajectory_file 
    """

    def __init__(self, settings_path: str, trajectory_file: str):
        try:
            self.traj = Trajectory(trajectory_file)
            self.settings = json.loads(Path(settings_path).read_text())
            self.elastic_properties = self._cubicConstantsFromTrajectory(ref=self.traj[0])
        except FileNotFoundError:
            raise FileNotFoundError(f"Trajectory file {trajectory_file} not found")

    def vizualize(self):
        view(self.traj)

    def computeCohesiveEnergy(self, discard_fraction: float = 0.2, return_unit: str = "J_per_atom") -> float:
        """
        Calculate the cohesive energy per atom.
        """

        calc = EMT()
        for j in range(len(self.traj[0])):
            atom = Atoms(self.traj[0].get_chemical_symbols()[j], positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=False)
            atom.calc = calc
            e_atom = atom.get_potential_energy()
        e_coh_list = []
        for i in range(len(self.traj)):
            e_bulk = self.traj[i].get_potential_energy()
            number = self.traj[i].get_global_number_of_atoms()
            e_coh_list.append(e_atom - e_bulk/number)

        e_coh_mean = np.mean(e_coh_list)
        logger.info(f"Cohesive energy: {e_coh_mean} eV")
        return e_coh_mean


    def computeLatticeConstant(self) -> float:
        """
        Calculates the lattice constant for fcc/bcc using atomic units internally.
        Returns the lattice constant in meters (SI).
        """
        # per-atom volume in Å^3
        V_A3_per_atom = self.traj[0].get_volume() / len(self.traj[0])

        struct = self.determineCrystalStructure()
        logger.info(f"Structure: {struct}")

        if struct == "fcc":
            a = (4.0 * V_A3_per_atom) ** (1.0 / 3.0)
        elif struct == "bcc":
            a = (2.0 * V_A3_per_atom) ** (1.0 / 3.0)
        else:
            # fallback: estimate from volume assuming simple cubic
            a = (1.0 * V_A3_per_atom) ** (1.0 / 3.0)

        logger.info(f"Lattice constant: {a} Å")
        return a

    def computeBulkModulus(self):
        """
        Computes bulk modulus from simulation values using atomic-unit processing.
        Returns bulk modulus in Pascals (SI). Only works if purposeful distortion is introduced in the simulation.
        """
        elastic_properties = self.elastic_properties
        K = elastic_properties['K']
        return float(K)

    def computeInternalPressure(self):
        """
        Compute internal pressure using atomic units internally.
        Instantaneous: P = (1/3V) [ 2 N E_kin + sum_i r_i · f_i ]
        where energies are in Hartree and lengths in Bohr; V in Bohr^3, P in Eh/a0^3.
        Return the time-average in SI Pascals.
        """
        internal_pressures_eVA3 = []
        for atoms in self.traj:
            e_kin_eV = atoms.get_kinetic_energy()
            N = len(atoms)
            V_A3 = atoms.get_volume()
            forces_eVA = atoms.get_forces()
            positions_A = atoms.get_positions()



            sum_rf_Eh = float(np.sum(forces_eVA * positions_A))

            P_eVA3 = (1.0 / (3.0 * V_A3)) * (2.0 * N * e_kin_eV + sum_rf_Eh)
            internal_pressures_eVA3.append(P_eVA3)

        avg_P = float(np.mean(internal_pressures_eVA3)) if internal_pressures_eVA3 else float('nan')
        logger.info(f"Average internal pressure = {avg_P} eV/Å3")
        return avg_P

    def computeMSD(self, time=-10, reference=0, return_SI=True):
        # TODO Shouldn't this be some mean value of many of the late snapshots in trajectory
        traj = self.traj
        if time is None:
            time = reference + 1

        N = traj[reference].get_global_number_of_atoms()
        r_0 = traj[reference].get_positions()  # Å

        r_n = traj[time].get_positions()  # Å

        msd = float(np.sum((r_n - r_0) ** 2) / N)
        logger.info(f"MSD : {msd} Å^2")
        if return_SI:
            msd_m2 = msd * (A_TO_M ** 2)
            logger.info(f"MSD (SI): {msd_m2} m^2")
            return msd_m2
        else:
            return msd

    def computeLindemannCriterion(self):
        """
        Returns True if Lindemann criterion is met, False otherwise.
        """

        cutoffs = natural_cutoffs(self.traj[0])
        neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
        neighbor_list.update(self.traj[0])
        min_dist = np.min([nearestNeighborDistanceForAtom(j, self.traj[0], neighbor_list) for j in
                           range(self.traj[0].get_global_number_of_atoms())])
        min_dist = float(min_dist)
        logger.info(f"First min dist-------------- {min_dist}")

        for atoms in self.traj:
            neighbor_list.update(atoms)
            nn_per_atom = np.min([nearestNeighborDistanceForAtom(i, atoms, neighbor_list) for i in
                                  range(atoms.get_global_number_of_atoms())])
            nn_per_atom = float(nn_per_atom)
            if min_dist is None:
                min_dist = nn_per_atom
            elif min_dist is not None and (nn_per_atom >= 0) and (nn_per_atom <= min_dist):
                min_dist = nn_per_atom

        logger.info(f"Overall nearest-neighbor distance [Å]: {min_dist}")
        # TODO ComputeMSD is only a snapshot, not mean. This is a problem I believe
        msd_m2 = self.computeMSD()  # SI m^2
        min_dist_m = min_dist * A_TO_M
        l = float(np.sqrt(msd_m2) / min_dist_m)
        logger.info(f"Lindemann index: {l}")
        if l >= 0.1:
            logger.info("Lindemann index is high.")
            return True
        else:
            logger.info("Lindemann index is low.")
            return False

    def computeSelfDiffusionCoefficient(self):  # Needs constant temperature, for current implementation
        timestep_list = []
        temp_list = [round(self.traj[i].get_temperature(), -1) for i in range(len(self.traj))]
        if self.settings["Temperature"] in temp_list:
            for i in range(len(self.traj)):
                if temp_list[i] == self.settings["Temperature"]:
                    timestep = i * self.settings["Timestep"] * self.settings["Interval"]
                    timestep_list.append([timestep, i])
        else:
            for i in range(len(self.traj)):
                timestep = i * self.settings["Timestep"] * self.settings["Interval"]
                timestep_list.append([timestep, i])

        # Compute MSD in atomic units
        msd0 = self.computeMSD(time=timestep_list[0][1], return_SI=False)
        msd1 = self.computeMSD(time=timestep_list[-1][1], return_SI=False)
        # Time difference
        t0_ps, t1_ps = timestep_list[0][0], timestep_list[-1][0]
        dt_fs = (t1_ps - t0_ps) * 1000

        slope = (msd1 - msd0) / dt_fs
        D = slope / 6.0
        logger.info(f"Self-Diffusion Coefficient: {D} A^2/fs")
        # Convert to SI m^2/s
        D_m2_per_s = D * (A_TO_M ** 2) * 1e-15
        logger.info(f"Self-Diffusion Coefficient: {D_m2_per_s} m^2/s")
        return float(D)

    def computeSpecificHeat(self):  # Requires NVT, might implement for NVE as well

        # total energy per frame in eV
        energy_eV = np.array(
            [atom_frame.get_potential_energy() + atom_frame.get_kinetic_energy() for atom_frame in self.traj])
        temperature = float(np.mean([atom_frame.get_temperature() for atom_frame in self.traj]))  # K

        e_mean = float(np.mean(energy_eV))
        e_2_mean = float(np.mean(energy_eV ** 2))
        prefactor = 1.0 / (kB * temperature ** 2)
        Cv_system = prefactor * (e_2_mean - e_mean ** 2)

        # Total mass in kg
        total_mass_amu = float(sum(self.traj[0].get_masses()))  # ASE masses in amu

        Cv_specific_eV_per_AmuK = Cv_system / total_mass_amu

        Cv_specific_J_per_kgK = Cv_specific_eV_per_AmuK * EV_TO_JOULE/ AMU_TO_KG
        logger.info(f"Specific heat capacity: {Cv_specific_J_per_kgK} J/(kg·K)")
        return float(Cv_specific_J_per_kgK)

    def computeDebyeTemperature(self):
        """
        Compute Debye temperature using atomic units internally.
        Returns Theta_D in Kelvin (SI).
        """

        if self.settings["Temperature"] < 1:
            # Low-temperature Debye from heat capacity; use a.u. for kB
            C_v = self.computeSpecificHeat()  # [J kg-1 K-1]
            temperature = float(np.mean([atom_frame.get_temperature() for atom_frame in self.traj]))
            N = self.traj[0].get_global_number_of_atoms()
            # Here we keep the classical constant form but ensure SI at the end
            debye = (234 * N * kB * EV_TO_JOULE * temperature ** 3 / C_v) ** (1 / 3)
            logger.info(f"Debye temperature: {debye} K")
            return float(debye)
        else:
            # Θ_D = (ħ/kB) (6π^2 n)^(1/3) vm
            out = self.elastic_properties  # SI Pa
            G = out['G']
            G_Pa = G * EV_PER_A3_TO_PA
            K = out['K']
            logger.info(f"Bulk modulus = {K * EV_PER_A3_TO_GPA} GPa")
            K_Pa = K * EV_PER_A3_TO_PA

            if K < 0 or G < 0:
                error_log = "Negative pressure during calulation of Debye temperature"
                logger.error(error_log)
                raise ValueError(error_log)

            if not (np.isfinite(G) and np.isfinite(K)):
                logger.info("Elastic constants not reliable from traj; Debye temperature cannot be computed.")
                return float('nan')


            V_A3 = np.mean([fr.get_volume() for fr in self.traj])
            if V_A3 is None:
                V_A3 = float(self.traj[0].get_volume())
            mass_u = float(sum(self.traj[0].get_masses()))
            rho = (mass_u / V_A3) * AMU_TO_KG * 1e30

            vT = math.sqrt(G_Pa / rho)
            vL = math.sqrt((K_Pa + 4.0 * G_Pa / 3.0) / rho)
            vm = ((1.0 / 3.0) * (1.0 / (vL ** 3) + 2.0 / (vT ** 3))) ** (-1.0 / 3.0)

            N = self.traj[0].get_global_number_of_atoms()
            n = (N / V_A3) * 1e30

            Theta_D = (hbar / kB) * vm * (6.0 * np.pi ** 2 * n) ** (1.0 / 3.0)
            logger.info(f"Debye temperature: {Theta_D} K")
            return float(Theta_D)

    def _cubicConstantsFromTrajectory(self, ref=None, tol_abs=1e-8, tol_rel=1e-3):
        """
        Estimate C11, C12, C44, K, G using only cell geometry and info['stress'] saved in traj.
        Internally uses atomic units (stress in Eh/a0^3), and returns SI Pascals.
        """
        if ref is None:
            ref = self.traj[0]

        # Collect arrays
        x_D, y_D = [], []  # for D = C11−C12 via (σa−σb) vs (εa−εb)
        x_S, y_S = [], []  # for S = C11+2C12 via σ_h vs tr(ε)
        x_xy, y_xy = [], []  # for C44 via σ_xy vs γ_xy
        x_xz, y_xz = [], []
        x_yz, y_yz = [], []

        for fr in self.traj:
            sig_eVA3 = _stressEVA3FromInfo(fr)
            if sig_eVA3 is None:
                continue
            eps = _smallStrainFromCells(ref, fr)

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

    def determineCrystalStructure(self):
        atoms = self.traj[0]
        [a, b, c, ang_bc, ang_ac, ang_ab] = atoms.get_cell_lengths_and_angles()
        ang_bc = round(ang_bc, 0)
        ang_ac = round(ang_ac, 0)
        ang_ab = round(ang_ab, 0)
        if (a == b and a == c and b == c) and (
                ang_bc == 90 and ang_ac == 90 and ang_ab == 90):  # Assumes primitive cell in POSCAR
            return "bcc"
        elif (a == b and a == c and b == c) and (
                ang_bc == 60 and ang_ac == 60 and ang_ab == 60):  # Assumes primitive cell in POSCAR
            return "fcc"

    def nearestNeightborsMean(self):
        """Calculate the mean distance of nearest neighbor in the structure for the last ten states of the simulation
        Loop structure: Last ten states -> Each atom -> neighbors to current atom
        """
        INF = 1e9
        NN_list = []
        for state in range(-10, -1):
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
        logger.info(f"Mean value of nearest neighbor : {NN_mean_distance}")
        return NN_mean_distance


def nearestNeighborDistanceForAtom(i, atoms, nl):
    indices, offsets = nl.get_neighbors(i)
    if len(indices) == 0:
        return float("nan")

    cell = atoms.cell.array  # (3, 3)
    translations = offsets @ cell  # (M, 3), periodic image translations
    rij = atoms.positions[indices] + translations - atoms.positions[i]  # (M, 3)
    dists = np.linalg.norm(rij, axis=1)

    # Remove exact or near-zero self entries if any slipped in
    positive = dists[dists > 1e-8]
    return float(np.min(positive)) if positive.size else float("nan")


def parseLogFile(log_path: str = "Outputs/output.log"):
    """
    Read data.log and return numeric columns.
    Keys: "Time[ps]", "Etot/N[eV]", "Epot/N[eV]", "Ekin/N[eV]", "T[K]".
    Repeated headers and bad lines are ignored.
    """
    text = Path(log_path).read_text(errors="ignore")
    lines = text.splitlines()
    headers = ["Time[ps]", "Etot/N[eV]", "Epot/N[eV]", "Ekin/N[eV]", "T[K]"]
    data = {h: [] for h in headers}
    for line in lines:
        line = line.strip()
        if not line or line.startswith("Time[ps]"):
            continue
        parts = re.split(r"\s+", line)
        if len(parts) < 5:
            continue
        try:
            t, etot, epot, ekin, T = (float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]),
                                      float(parts[4]))
        except ValueError:
            continue
        # Store
        data[headers[0]].append(t)
        data[headers[1]].append(etot)
        data[headers[2]].append(epot)
        data[headers[3]].append(ekin)
        data[headers[4]].append(T)
    return data


def plot(msd_list, dt=1.0, title="Mean Squared Displacement", save_as=None):
    import matplotlib.pyplot as plt

    """
    Plot MSD vs time.

    Parameters
    ----------
    msd_list : list or np.ndarray
        MSD values per frame.
    dt : float
        Timestep size in fs (default 1 fs).
    title : str
        Plot title.
    save_as : str or None
        Filename to save the plot. If None, plot is shown interactively.
    """
    msd_array = np.array(msd_list)
    time = np.arange(1, len(msd_array) + 1) * dt  # start from frame 1

    plt.figure(figsize=(6, 4))
    plt.plot(time, msd_array, marker='o', linestyle='-', color='tab:blue')
    plt.xlabel("Time [fs]")
    plt.ylabel("MSD [Å²]")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()

    if save_as:
        plt.savefig(save_as, dpi=300)
        plt.close()
    else:
        plt.show()


def _smallStrainFromCells(ref, cur):
    C_ref = ref.cell.array
    C_cur = cur.cell.array
    A_ref = C_ref.T
    #logger.info(f"A_ref: {A_ref} \n")
    A_cur = C_cur.T
    #logger.info(f"A_cur: {A_cur} \n")
    #logger.info(f"INV : {np.linalg.inv(A_ref)} \n")
    F = A_cur @ np.linalg.inv(A_ref)
    #logger.info(f"F: {F} \n")
    eps = 0.5 * (F + F.T) - np.eye(3)
    return eps


def _stressEVA3FromInfo(cur):
    sig_v = cur.info["stress eV/A3"]
    if sig_v is None:
        return None
    return np.asarray(sig_v, dtype=float)  # eV/A3
