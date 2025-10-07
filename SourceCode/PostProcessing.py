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

# Common conversions
EV_PER_A3_TO_GPA = 160.21766208
EV_PER_A3_TO_PA = 160.21766208e9  # Pa per (eV/Å^3)
AMU_TO_KG = 1.66053906660e-27
EV_TO_JOULE = 1.602176634e-19
AVOGRADO = 6.02214076e23
A_TO_M = 1e-10

# Atomic unit constants
BOHR_ANG = Bohr                 # 1 a0 in Å
BOHR_M = BOHR_ANG * A_TO_M      # 1 a0 in m
HARTREE_eV = Hartree            # 1 Eh in eV
HARTREE_J = HARTREE_eV * EV_TO_JOULE  # 1 Eh in J
AU_TIME_S = physical_constants['atomic unit of time'][0]  # s
AU_VEL_MS = BOHR_M / AU_TIME_S  # m/s
AU_PRESSURE_PA = HARTREE_J / (BOHR_M**3)  # Pa
AMU_TO_ME = physical_constants['atomic mass constant energy equivalent in MeV'][0]  # placeholder, replace below
# Better: ratio amu to electron mass
AMU_TO_ME = physical_constants['atomic mass constant'][0] / physical_constants['electron mass'][0]

# Helper conversions to/from atomic units
# Length
ang_to_bohr = lambda x: np.asarray(x, dtype=float) / BOHR_ANG
bohr_to_m = lambda x: np.asarray(x, dtype=float) * BOHR_M
# Energy
ev_to_hartree = lambda x: np.asarray(x, dtype=float) / HARTREE_eV
hartree_to_j = lambda x: np.asarray(x, dtype=float) * HARTREE_J
# Volume
A3_to_bohr3 = lambda x: np.asarray(x, dtype=float) / (BOHR_ANG**3)
# Force
eV_per_A_to_Eh_per_bohr = lambda x: (np.asarray(x, dtype=float) / HARTREE_eV) * (1.0 / BOHR_ANG)
# Pressure
eV_per_A3_to_auP = lambda x: (np.asarray(x, dtype=float) / HARTREE_eV) * (BOHR_ANG**3)
auP_to_Pa = lambda x: np.asarray(x, dtype=float) * AU_PRESSURE_PA

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
            self.elastic_properties = self._cubic_constants_from_traj(ref=self.traj[0])
        except FileNotFoundError:
            raise FileNotFoundError(f"Trajectory file {trajectory_file} not found")

    def vizualize(self):
        view(self.traj)

    def computeCohesiveEnergy(self, discard_fraction: float = 0.2, return_unit: str = "J_per_atom") -> float:
        """
        Estimate cohesive energy per atom using only simulation data (no external calculator).
        Principle (harmonic solid): U(T) = U0 + K(T) ⇒ U0 = U(T) − K(T).
        If the potential zero is at infinite separation (true for common MD potentials),
        E_coh = −U0/N = K/N − U/N.

        Internals use atomic units (Hartree); returns SI Joules per atom by default.

        Parameters
        ----------
        discard_fraction : float
            Fraction of initial frames to discard as transient (default 0.2).
        return_unit : str
            "J_per_atom" (default) or "eV_per_atom".
        """

        pot_list = []
        kin_list = []
        N = int(self.traj[0].info['number_of_atoms'])
        for snapshot in self.traj:
            try:
                pot = float(snapshot.info['potential_energy eV']) / N
                kin = float(snapshot.info['kinetic_energy eV']) / N
                pot_list.append(pot)
                kin_list.append(kin)
            except Exception:
                logger.error("Could not acquire energy per atom for the cohesive energy")
                raise RuntimeError("Could not acquire energy per atom for the cohesive energy")
        Epot_per_atom_eV = np.asarray(pot_list, dtype=float)
        Ekin_per_atom_eV = np.asarray(kin_list, dtype=float)

        # Discard an initial fraction to avoid transients
        start = int(round(discard_fraction * Epot_per_atom_eV.size)) if Epot_per_atom_eV.size > 5 else 0
        EpotN_eV = Epot_per_atom_eV[start:]
        EkinN_eV = Ekin_per_atom_eV[start:]

        # Convert to atomic units and compute per-frame cohesive energy estimates
        EpotN_Eh = ev_to_hartree(EpotN_eV)
        EkinN_Eh = ev_to_hartree(EkinN_eV)
        Ecoh_Eh_per_atom_series = EkinN_Eh - EpotN_Eh  # = -U0/N in Eh

        # Average and convert at return
        Ecoh_Eh_per_atom = float(np.mean(Ecoh_Eh_per_atom_series))
        Ecoh_eV_per_atom = float(Ecoh_Eh_per_atom * HARTREE_eV)
        Ecoh_J_per_atom = float(hartree_to_j(Ecoh_Eh_per_atom))

        logger.info(f"Cohesive energy (from simulation): {Ecoh_eV_per_atom} eV/atom ≈ {Ecoh_J_per_atom} J/atom")

        if return_unit == "eV_per_atom":
            return Ecoh_eV_per_atom
        else:
            return Ecoh_J_per_atom

    def computeLatticeConstant(self) -> float:
        """
        Calculates the lattice constant for fcc/bcc using atomic units internally.
        Returns the lattice constant in meters (SI).
        """
        # per-atom volume in Å^3
        V_A3_per_atom = self.traj[0].info['volume A3'] / self.traj[0].info['number_of_atoms']
        V_bohr3 = A3_to_bohr3(V_A3_per_atom)

        struct = self.determineCrystalStructure()
        logger.info(f"Structure: {struct}")

        if struct == "fcc":
            a_bohr = (4.0 * V_bohr3) ** (1.0 / 3.0)
        elif struct == "bcc":
            a_bohr = (2.0 * V_bohr3) ** (1.0 / 3.0)
        else:
            # fallback: estimate from volume assuming simple cubic
            a_bohr = (1.0 * V_bohr3) ** (1.0 / 3.0)

        a_m = float(bohr_to_m(a_bohr))
        logger.info(f"Lattice constant: {a_m} m (SI)")
        return a_m

    def computeBulkModulus(self):
        """
        Computes bulk modulus from simulation values using atomic-unit processing.
        Returns bulk modulus in Pascals (SI). Only works if purposeful distortion is introduced in the simulation.
        """
        elastic_properties = self.elastic_properties
        logger.info(f"Elasticity properties (SI): {elastic_properties}")
        K_Pa = elastic_properties.get('K_Pa', float('nan'))
        logger.info(f"Bulk modulus: {K_Pa} Pa ({K_Pa*1e-9 if np.isfinite(K_Pa) else float('nan')} GPa)")
        return float(K_Pa)

    def computeInternalPressure(self):
        """
        Compute internal pressure using atomic units internally.
        Instantaneous: P = (1/3V) [ 2 N E_kin + sum_i r_i · f_i ]
        where energies are in Hartree and lengths in Bohr; V in Bohr^3, P in Eh/a0^3.
        Return the time-average in SI Pascals.
        """
        internal_pressures_Pa = []
        for atoms in self.traj:
            e_kin_eV = atoms.info['kinetic_energy eV']
            N = atoms.info['number_of_atoms']
            V_A3 = atoms.info['volume A3']
            forces_eVA = atoms.info['forces eV/A']
            positions_A = atoms.info['positions']

            # Convert to atomic units
            Ekin_Eh = ev_to_hartree(e_kin_eV)
            V_bohr3 = A3_to_bohr3(V_A3)
            F_Eh_per_bohr = eV_per_A_to_Eh_per_bohr(forces_eVA)
            R_bohr = ang_to_bohr(positions_A)
            sum_rf_Eh = float(np.sum(F_Eh_per_bohr * R_bohr))

            P_au = (1.0 / (3.0 * V_bohr3)) * (2.0 * N * Ekin_Eh + sum_rf_Eh)
            P_Pa = float(auP_to_Pa(P_au))
            internal_pressures_Pa.append(P_Pa)

        avg_Pa = float(np.mean(internal_pressures_Pa)) if internal_pressures_Pa else float('nan')
        logger.info(f"Average internal pressure = {avg_Pa} Pa")
        return avg_Pa

    def computeMSD(self, time=-10, reference=0, return_SI=True):
        # TODO Shouldn't this be some mean value of many of the late snapshots in trajectory
        traj = self.traj
        if time is None:
            time = reference + 1

        N = traj[reference].get_global_number_of_atoms()
        r_0 = traj[reference].get_positions()  # Å

        r_n = traj[time].get_positions()  # Å
        # convert to atomic units (Bohr), do calculation in a.u.
        r0_bohr = ang_to_bohr(r_0)
        rn_bohr = ang_to_bohr(r_n)
        msd_bohr2 = float(np.sum((rn_bohr - r0_bohr) ** 2) / N)
        logger.info(f"MSD (a.u.): {msd_bohr2} a0^2")
        if return_SI:
            msd_m2 = msd_bohr2 * (BOHR_M ** 2)
            logger.info(f"MSD (SI): {msd_m2} m^2")
            return msd_m2
        else:
            return msd_bohr2


    def computeLindemannCriterion(self):
        """
        Returns True if Lindemann criterion is met, False otherwise.
        """

        cutoffs = natural_cutoffs(self.traj[0])
        neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
        neighbor_list.update(self.traj[0])
        min_dist = np.min([nearest_neighbor_distance_for_atom(j, self.traj[0], neighbor_list) for j in range(self.traj[0].get_global_number_of_atoms())])
        min_dist = float(min_dist)
        logger.info(f"First min dist-------------- {min_dist}")

        for atoms in self.traj:
            neighbor_list.update(atoms)
            nn_per_atom = np.min([nearest_neighbor_distance_for_atom(i, atoms, neighbor_list) for i in range(atoms.get_global_number_of_atoms())])
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
        temp_list = [round(self.traj[i].info["temperature"], -1) for i in range(len(self.traj))]
        if self.settings["Temperature"] in temp_list:
            for i in range(len(self.traj)):
                if temp_list[i] == self.settings["Temperature"]:
                    timestep = i * self.settings["Timestep"] * self.settings["Interval"]
                    timestep_list.append([timestep, i])
        else:
            for i in range(len(self.traj)):
                timestep = i * self.settings["Timestep"] * self.settings["Interval"]
                timestep_list.append([timestep, i])

        logger.info(f"{timestep_list[0][0]} -------- {timestep_list[-1][0]}")

        # Compute MSD in atomic units (Bohr^2)
        msd0_bohr2 = self.computeMSD(time=timestep_list[0][1], return_SI=False)
        msd1_bohr2 = self.computeMSD(time=timestep_list[-1][1], return_SI=False)
        # Time difference in atomic units
        t0_ps, t1_ps = timestep_list[0][0], timestep_list[-1][0]
        dt_s = (t1_ps - t0_ps) * 1e-12
        dt_au = dt_s / AU_TIME_S

        slope_bohr2_per_au = (msd1_bohr2 - msd0_bohr2) / dt_au
        D_bohr2_per_au = slope_bohr2_per_au / 6.0
        # Convert to SI m^2/s
        D_m2_per_s = D_bohr2_per_au * (BOHR_M**2) / AU_TIME_S
        logger.info(f"Self-Diffusion Coefficient: {D_m2_per_s} m^2/s")
        return float(D_m2_per_s)

    def computeSpecificHeat(self):  # Requires NVT, might implement for NVE as well

        # total energy per frame in eV; convert to Hartree for a.u. computation
        energy_eV = np.array([atom_frame.get_potential_energy() + atom_frame.get_kinetic_energy() for atom_frame in self.traj])
        energy_Eh = ev_to_hartree(energy_eV)
        temperature = float(np.mean([atom_frame.get_temperature() for atom_frame in self.traj]))  # K

        e_mean = float(np.mean(energy_Eh))    # [Eh]
        e_2_mean = float(np.mean(energy_Eh ** 2))     # [Eh^2]
        # Boltzmann constant in Hartree/K
        kB_Eh_per_K = kB / HARTREE_eV
        prefactor = 1.0 / (kB_Eh_per_K * temperature**2)   # [1 / (K Eh)]
        Cv_system_Eh_per_K = prefactor * (e_2_mean - e_mean**2)     # [Eh/K]
        Cv_system_J_per_K = float(hartree_to_j(Cv_system_Eh_per_K))

        # Total mass in kg
        total_mass_amu = float(sum(self.traj[0].get_masses()))  # ASE masses in amu
        mass_kg = total_mass_amu * AMU_TO_KG

        Cv_specific_J_per_kgK = Cv_system_J_per_K / mass_kg
        logger.info(f"Specific heat capacity: {Cv_specific_J_per_kgK} J/(kg·K)")
        return float(Cv_specific_J_per_kgK)

    def computeDebyeTemperature(self):
        """
        Compute Debye temperature using atomic units internally.
        Returns Theta_D in Kelvin (SI).
        """

        if self.settings["Temperature"] < 1:
            # Low-temperature Debye from heat capacity; use a.u. for kB
            C_v = self.computeSpecificHeat()    # [J kg-1 K-1]
            temperature = float(np.mean([atom_frame.get_temperature() for atom_frame in self.traj]))
            N = self.traj[0].get_global_number_of_atoms()
            # Here we keep the classical constant form but ensure SI at the end
            debye = (234 * N * kB * EV_TO_JOULE * temperature**3 / C_v) ** (1/3)
            logger.info(f"Debye temperature: {debye} K")
            return float(debye)
        else:
            # Θ_D = (ħ/kB) (6π^2 n)^(1/3) vm
            out = self.elastic_properties  # SI Pa
            G_Pa = out['G_Pa']
            K_Pa = out['K_Pa'] if np.isfinite(out['K_Pa']) else (3.0 * out['C11_Pa'] - 2.0 * out['C44_Pa']) / 3.0

            if not (np.isfinite(G_Pa) and np.isfinite(K_Pa)):
                logger.info("Elastic constants not reliable from traj; Debye temperature cannot be computed.")
                return float('nan')

            # Convert elastic constants to atomic units of pressure
            G_au = G_Pa / AU_PRESSURE_PA
            K_au = K_Pa / AU_PRESSURE_PA

            # Density in atomic units (electron masses per Bohr^3)
            V_A3 = np.mean([fr.info['volume A3'] for fr in self.traj])
            if V_A3 is None:
                V_A3 = float(self.traj[0].get_volume())
            V_bohr3 = A3_to_bohr3(V_A3)
            mass_u = float(sum(self.traj[0].get_masses()))
            mass_me = mass_u * AMU_TO_ME
            rho_au = mass_me / V_bohr3  # m_e / a0^3

            # Sound velocities in atomic units (a0 / au_time)
            vT_au = math.sqrt(G_au / rho_au)
            vL_au = math.sqrt((K_au + 4.0 * G_au / 3.0) / rho_au)
            vm_au = ((1.0 / 3.0) * (1.0 / (vL_au ** 3) + 2.0 / (vT_au ** 3))) ** (-1.0 / 3.0)

            # Atom number density in a.u. (per Bohr^3)
            N = self.traj[0].get_global_number_of_atoms()
            n_au = N / V_bohr3

            # In atomic units, ħ = 1; kB in Hartree/K
            kB_Eh_per_K = kB / HARTREE_eV
            Theta_D = (1.0 / kB_Eh_per_K) * vm_au * (6.0 * pi ** 2 * n_au) ** (1.0 / 3.0)
            logger.info(f"Debye temperature: {Theta_D} K")
            return float(Theta_D)


    def _cubic_constants_from_traj(self, ref=None, tol_abs=1e-8, tol_rel=1e-3):
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
            sig_eVA3 = _stress_eVA3_from_info(fr)
            if sig_eVA3 is None:
                continue
            # convert stress to atomic units (Eh/a0^3)
            sig_au = eV_per_A3_to_auP(sig_eVA3)
            eps = _small_strain_from_cells(ref, fr)

            e_xx, e_yy, e_zz = eps[0, 0], eps[1, 1], eps[2, 2]
            e_xy, e_xz, e_yz = eps[0, 1], eps[0, 2], eps[1, 2]
            g_xy, g_xz, g_yz = 2 * e_xy, 2 * e_xz, 2 * e_yz

            s_xx, s_yy, s_zz = sig_au[0], sig_au[1], sig_au[2]
            s_yz, s_xz, s_xy = sig_au[3], sig_au[4], sig_au[5]

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

        def slope_filtered(x, y):
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

        D_au = slope_filtered(x_D, y_D)
        S_3_au = slope_filtered(x_S, y_S)
        S_3_au = S_3_au * 3.0 if np.isfinite(S_3_au) else float('nan')

        C44_xy_au = slope_filtered(x_xy, y_xy)
        C44_xz_au = slope_filtered(x_xz, y_xz)
        C44_yz_au = slope_filtered(x_yz, y_yz)

        C44s_au = [v for v in [C44_xy_au, C44_xz_au, C44_yz_au] if np.isfinite(v)]
        C44_au = float(np.mean(C44s_au)) if C44s_au else float('nan')

        if np.isfinite(S_3_au) and np.isfinite(D_au):
            C11_au = (S_3_au + 2.0 * D_au) / 3.0
            C12_au = (S_3_au - D_au) / 3.0
            K_au = S_3_au / 3.0
        else:
            C11_au = C12_au = K_au = float('nan')

        # VRH shear modulus
        if np.isfinite(D_au) and np.isfinite(C44_au):
            G_V_au = (D_au + 3.0 * C44_au) / 5.0
            denom = 4.0 * C44_au + 3.0 * D_au
            if denom != 0:
                G_R_au = 5.0 * D_au * C44_au / denom
                G_au = 0.5 * (G_V_au + G_R_au) if np.isfinite(G_R_au) else G_V_au
            else:
                G_au = G_V_au
        else:
            G_au = C44_au

        # Convert to SI Pa
        C11_Pa = float(auP_to_Pa(C11_au)) if np.isfinite(C11_au) else float('nan')
        C12_Pa = float(auP_to_Pa(C12_au)) if np.isfinite(C12_au) else float('nan')
        C44_Pa = float(auP_to_Pa(C44_au)) if np.isfinite(C44_au) else float('nan')
        K_Pa = float(auP_to_Pa(K_au)) if np.isfinite(K_au) else float('nan')
        G_Pa = float(auP_to_Pa(G_au)) if np.isfinite(G_au) else float('nan')

        return {
            'C11_Pa': C11_Pa, 'C12_Pa': C12_Pa, 'C44_Pa': C44_Pa, 'K_Pa': K_Pa, 'G_Pa': G_Pa,
            'C44_xy_Pa': float(auP_to_Pa(C44_xy_au)) if np.isfinite(C44_xy_au) else float('nan'),
            'C44_xz_Pa': float(auP_to_Pa(C44_xz_au)) if np.isfinite(C44_xz_au) else float('nan'),
            'C44_yz_Pa': float(auP_to_Pa(C44_yz_au)) if np.isfinite(C44_yz_au) else float('nan')
        }

    def determineCrystalStructure(self):
        atoms = self.traj[0]
        [a, b, c, ang_bc, ang_ac, ang_ab] = atoms.get_cell_lengths_and_angles()
        ang_bc = round(ang_bc, 0)
        ang_ac = round(ang_ac, 0)
        ang_ab = round(ang_ab, 0)
        if (a == b and a == c and b == c) and (ang_bc == 90 and ang_ac == 90 and ang_ab == 90): # Assumes primitive cell in POSCAR
            return "bcc"
        elif (a == b and a == c and b == c) and (ang_bc == 60 and ang_ac == 60 and ang_ab == 60):   # Assumes primitive cell in POSCAR
            return "fcc"

    def NearestNeightborsMean(self, trajectory_index: int = 0):
        """Calculate the mean distance of nearest neighbor in the strcuture"""
        atoms = self.traj[trajectory_index]
        for i in range(atoms.get_global_number_of_atoms()):
            for j in range(i, atoms.get_global_number_of_atoms()):
                cutoff = natural_cutoffs(atoms)
                neighbor_list = NeighborList(cutoff, self_interaction=False, bothways=True)
                neighbor_list.update(atoms)



def nearest_neighbor_distance_for_atom(i, atoms, nl):
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


def ParseLogFile(log_path: str = "Outputs/output.log"):
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
            t, etot, epot, ekin, T = (float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]))
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



def _small_strain_from_cells(ref, cur):
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


def _stress_eVA3_from_info(cur):
    sig_v = cur.info['stress eV/A3']
    if sig_v is None:
        return None
    return np.asarray(sig_v, dtype=float) #eV/A3