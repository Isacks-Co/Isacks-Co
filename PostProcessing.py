from random import randint

from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.io import read, write
from ase import Atoms
import re
from pathlib import Path
import numpy as np
import json
from ase.eos import EquationOfState
from ase.calculators.emt import EMT
from ase.calculators.eam import EAM
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list, NeighborList
import json
import mpmath as mp
from ase.md.analysis import DiffusionCoefficient
from ase.units import Bohr, Rydberg, kJ, kB, fs as fs_conversion, Hartree, mol, kcal, Ang
from ase.build import bulk
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii, cohesive_energies
import MDAnalysis as mda
from MDAnalysis.analysis import msd as mda_msd
from MDAnalysis.tests.datafiles import RANDOM_WALK_TOPO, RANDOM_WALK
import matplotlib.pyplot as plt
from scipy.stats import linregress
import PreProcessing
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from scipy.constants import pi, hbar, k as kB
from phonopy.interface.calculator import read_crystal_structure
from math import sqrt

EV_PER_A3_TO_GPA = 160.21766208
ATOMIC_MASS_IN_KG = 1.66053906892e-27

class PostProcessing:
    """
        Uses the log files generates from the Molecular Dynamics class to vizualize and analyze the data
        input trajectory file as traj_file 
    """
    def __init__(self, pp, traj_file = "output.traj", settings_path = "settings.json", log_path = "output.log"):
        try:
            self.traj = Trajectory(traj_file)
            self.settings = json.loads(Path(settings_path).read_text())
            self.data_log = ParseLogFile(log_path)
            self.ideal_atoms = pp.readAtomicStructure('POSCAR').repeat( tuple([self.settings["Supercells"]+1]*3))
        except FileNotFoundError:
            raise FileNotFoundError(f"Trajectory file {traj_file} not found")
        
    def vizualize(self):
        view(self.traj)

    def computeCohesiveEnergy(self, poscar_path: str = "POSCAR", settings_path: str = "settings.json") -> float:
        """
        Compute cohesive energy per atom for the crystal in POSCAR.
        Definition: E_coh = (sum_i n_i E_i^atom − E_crystal_total) / N
        Note: Uses a fixed EMT calculator; does not read settings or apply fallbacks.
        """

        # Read structure
        atoms = read(poscar_path)

        # Taget från databas
        atomic_number = (atoms.get_atomic_numbers())[0]
        energy = cohesive_energies.cohesive_energy_kittel2005[atomic_number]
        print("E_coh _--------------------------------------------------= ", energy)
        return energy

        # Manuell uträkning
        """
        # Crystal total energy
        atoms_cryst = atoms.copy()
        atoms_cryst.calc = EMT()
        E_cryst = float(atoms_cryst.get_potential_energy())

        # Isolated atom reference energies per species
        symbols = atoms_cryst.get_chemical_symbols()
        unique = {}
        for s in set(symbols):
            one = Atoms(symbols=s, positions=[[0.0, 0.0, 0.0]], cell=[20.0, 20.0, 20.0], pbc=False)
            one.calc = EMT()
            unique[s] = float(one.get_potential_energy())

        # Sum over composition
        N = len(symbols)
        sum_iso = sum(unique[s] for s in symbols)

        # Cohesive energy per atom
        E_coh = (sum_iso - E_cryst) / N
        print("E_coh = ", E_coh)
        return float(E_coh)
        """

    def computeLatticeConstant(self, poscar_path: str = "POSCAR") -> float:
        """
        Calculates the lattice constant of a FCC crystal.
        """
        atoms = read(poscar_path)
        # 3x3 lattice matrix (rows are the lattice vectors in Cartesian Å)
        cell = atoms.cell.array  # numpy.ndarray of shape (3, 3)

        a1, a2, a3 = cell[0], cell[1], cell[2]
        # compute volume
        V = abs(np.dot(a1, np.cross(a2, a3)))

        struct = self.determineCrystalStructure(poscar_path)
        print("Structure: ", struct)

        if struct == "fcc":
            # If these are primitive fcc vectors:
            a_from_V = (4.0 * V) ** (1.0 / 3.0)
            print("a (from primitive-cell volume) = {:.6f} Å".format(a_from_V))
            return float(a_from_V)

        elif struct == "bcc":
            # If these are primitive bcc vectors:
            a_from_V = (2.0 * V) ** (1.0 / 3.0)
            print("a (from primitive-cell volume) = {:.6f} Å".format(a_from_V))
            return float(a_from_V)

        else:
            # If these were conventional cubic vectors (orthogonal, equal length):
            a_conv = np.linalg.norm(a1)  # would equal the conventional a
            print("a (if conventional cell) = {:.6f} Å".format(a_conv))
            #Added the section immediately above, because I don't remember much of solid state physics, and want to be sure I'm doing this right.
            return float(a_conv)

    def computeBulkModulus(self, poscar_path: str = "POSCAR", scales = np.linspace(0.96, 1.04, 7)):
        """
        Computes bulk modulus from equation of states, by rescaling.
        """
        atoms0 = read(poscar_path)

        vols, energies = [], []
        for s in scales:
            a = atoms0.copy()
            a.set_cell(atoms0.get_cell() * s, scale_atoms=True)
            a.calc = EMT()
            E = a.get_potential_energy()
            V = a.get_volume()
            energies.append(E)
            vols.append(V)

        eos = EquationOfState(vols, energies, eos='birchmurnaghan')
        v0, e0, B0 = eos.fit()  # B0 in eV/Å^3
        bulk_modulus = B0 * EV_PER_A3_TO_GPA
        #print("Bulk modulus = {:.6f} GPa".format(bulk_modulus))
        return float(bulk_modulus)

    def computeInternalPressure(self):
        """
        For NVT ensemble, computes internal pressure.
        Instantaneous: P(t)  = (1/3V)[2NkT(t) + SUM_i(r_i*f_i)]
        Average:       P = (1/M)SUM_i(P(nΔt))
        """
        traj = self.traj
        internal_pressure = []
        for atoms in traj:
            e_kin = atoms.get_kinetic_energy()
            number_of_atoms = atoms.get_global_number_of_atoms()
            volume = atoms.get_volume()
            internal_pressure.append(((1/(3*volume))*((2*number_of_atoms*e_kin) + np.sum(atoms.get_forces()*atoms.get_positions())))*EV_PER_A3_TO_GPA)
        internal_pressure = np.array(internal_pressure)
        average_internal_pressure = np.average(internal_pressure)
        print("Average internal pressure = ", average_internal_pressure, "GPa")
        return float(average_internal_pressure)

    def computeMSD(self, time = -10 , reference = 0, log_path = "output.log", settings_path = "settings.json"):
        traj = self.traj
        data_log = ParseLogFile(log_path)
        if time is None:
            time = reference + 1

        if reference is None:
            N = self.ideal_atoms.get_global_number_of_atoms()
            r_0 = self.ideal_atoms.get_positions()
        else:
            N = traj[reference].get_global_number_of_atoms()
            r_0 = traj[reference].get_positions()

        #print("r_0:---------------------- ", r_0)
        r_n = traj[time].get_positions()
        #print("r_n:---------------------- ", r_n)
        msd = np.sum((r_n - r_0)**2) / N        #Å^2
        print("MSD:---------------------- ", msd)
        return msd


    def computeLindemannCriterion(self, poscar_path: str = "POSCAR", settings_path = "settings.json"):
        """
        Returns True if Lindemann criterion is met, False otherwise.
        """
        cutoffs = natural_cutoffs(self.ideal_atoms)
        nl_initial = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl_initial.update(self.ideal_atoms)
        min_dist = np.min([nearest_neighbor_distance_for_atom(j, self.ideal_atoms, nl_initial) for j in range(self.ideal_atoms.get_global_number_of_atoms())])
        min_dist = float(min_dist)
        print("First min dist-------------- ", min_dist)

        traj = self.traj
        nl = build_neighbor_list(traj[0])
        for atoms in traj:
            nl.update(atoms)
            nn_per_atom = np.min([nearest_neighbor_distance_for_atom(i, atoms, nl) for i in range(atoms.get_global_number_of_atoms())])
            nn_per_atom = float(nn_per_atom)
            if min_dist is None:
                min_dist = nn_per_atom
            elif min_dist is not None and (nn_per_atom >= 0) and (nn_per_atom <= min_dist):
                min_dist = nn_per_atom

        print("Overall nearest-neighbor distance [Å]: ", min_dist)
        msd = self.computeMSD()
        l = np.sqrt(msd) / min_dist
        print("Lindemann index: ", l)
        if l >= 0.1:
            print("Lindemann index is high.")
            return True
        else:
            print("Lindemann index is low.")
            return False

    def computeSelfDiffusionCoefficient(self, settings_path="settings.json", log_path="output.log"):

        traj = self.traj
        data_log = ParseLogFile(log_path)
        try:
            settings = json.loads(Path(settings_path).read_text())
            timestep_fs = float(settings["Timestep"])  # fs per MD step
            interval = int(settings["Interval"])  # steps between saved frames
        except Exception as e:
            raise RuntimeError(f"Failed to read timing info from {settings_path}: {e}")

        timestep_list = []
        temp_list = [round(data_log["T[K]"][i], -1) for i in range(len(data_log["T[K]"]))]
        print("Temperature list: ", temp_list)
        if settings["Temperature"] in temp_list:
            for i in range(len(data_log["Time[ps]"])):
                if temp_list[i] == settings["Temperature"]:
                    timestep = data_log["Time[ps]"][i]
                    timestep_list.append([timestep, i])
        else:
            for i in range(len(data_log["Time[ps]"])):
                timestep = data_log["Time[ps]"][i]
                timestep_list.append([timestep, i])

        print("Timestep list: ", timestep_list)
        print(timestep_list[0][0], "--------", timestep_list[-1][0])


        #halfway = int(round((len(timestep_list) - 1) / 2, 0))
        msd_list = [self.computeMSD(time=timestep_list[0][1]),
                    self.computeMSD(time=timestep_list[-1][1]),]
        slope = (msd_list[-1] - msd_list[0])/ (timestep_list[-1][0] - timestep_list[0][0])
        slope_m2_per_s = slope * 1e-8 / 6
        print("Self-Diffusion Coefficient------------ [Å^2/s^-1]: ", slope_m2_per_s)
        return slope_m2_per_s


    def computeDebyeTemperature(self, log_path: str = "output.log", settings_path: str = "settings.json", poscar_path = "POSCAR" , flags = []):
        """
        atoms: ASE Atoms with an attached calculator that yields forces.
        supercell: supercell size for force constants.
        disp: displacement amplitude [Å].
        qmax_rlu: max fractional q (in reciprocal lattice units) for linear slope fit.

        Returns: dict with vT, vL, vm [m/s] and ThetaD [K].
        """

        traj = self.traj
        data_log = ParseLogFile(log_path)
        settings = json.loads(Path(settings_path).read_text())

        """
        # Θ_D = [12 * π⁴ * Nk_B * T³/ (5*C_v)] ^ (1/3)            only works when T<<Θ_D ~ <<1/C_v.     Note: 12 * π⁴ * Nk_B / 5 ≈ 234
        C_v = 385    # [J kg-1 K-1]     To be replaced with the return value from the C_v function, for the material
        D3 = []
        for i in range(len(data_log["T[K]"])):
            D3.append(234  * (data_log["T[K]"][i])**3 / C_v)
        print("Debye temperature ^3: ", D3)
        D3_mean = np.average(D3)
        print("Average Debye temperature ^3: ", D3_mean)
        Θ_D = round(D3_mean**(1/3), 0)
        print("Θ_D = ", Θ_D)
        """

        # Θ_D = (ħ/kB) (6π^2 n)^(1/3) vm
        debye = []
        for i in range(len(traj)):
            a = traj[i]
            #a.calc = EAM(potential="Cu_u3.eam.alloy")
            a.calc = EMT()
            V = a.get_volume() * 1e-30
            N = a.get_global_number_of_atoms()
            out = self.computeShearModulus_from_elastic_cubic(a)
            #print(out)  # G in GPa

            shear_modulus = out["G"] * 1e9
            mass = sum(a.get_masses())
            density = (mass / V) * ATOMIC_MASS_IN_KG
            vt = sqrt(shear_modulus / density)

            bulk_modulus = self.computeBulkModulus() * 1e9
            vl = sqrt((bulk_modulus + 4*shear_modulus/3)/density)

            vm = ((3/((1/(vl**3)) + (2/(vt**3))))**(1/3))

            n = N / V
            D = (hbar/kB) * vm * (6 * (pi**2) * n)**(1/3)
            #print("Debye temperature: ", D)
            debye.append(D)

        debye_avg = np.average(debye)
        print("Average Debye temperature: ", debye_avg)
        return debye_avg


    def computeShearModulus_from_elastic_cubic(self, atoms, eps_list=(2e-3, 3e-3, 5e-3)):
        """
        Compute cubic elastic constants (C11, C12, C44), bulk modulus K and shear modulus G (VRH).
        Intended for fcc (or any cubic) cells. Requires a calculator that provides stress.

        Returns dict with C11, C12, C44, K, G [GPa].
        """
        C0 = atoms.cell.array

        def stress_of(cell):
            a = atoms.copy()
            a.calc = atoms.calc  # ensure calculator on copy
            a.set_cell(cell, scale_atoms=True)
            s = a.get_stress(voigt=True) * EV_PER_A3_TO_GPA  # GPa
            return s

        eps_list = np.atleast_1d(eps_list)

        # 1) Orthorhombic (volume-conserving to 1st order) → C11 - C12
        #    Using σ1-σ2 ≈ 2 (C11 - C12) ε
        sig_diffs = []
        for eps in eps_list:
            Bp = np.diag([1 + eps, 1 - eps, 1 / (1 - eps ** 2)])
            Bm = np.diag([1 - eps, 1 + eps, 1 / (1 - eps ** 2)])
            sp = stress_of(C0 @ Bp)
            sm = stress_of(C0 @ Bm)
            sig = (sp - sm) / 2.0
            sig_diffs.append(0.5 * (sig[0] - sig[1]))  # ≈ (C11 - C12) * 2ε
        # Linear fit through origin: sig_diff ≈ 2 (C11 - C12) ε
        coeff = np.polyfit(eps_list, sig_diffs, 1)[0]
        C11_minus_C12 = 0.5 * coeff  # GPa

        # 2) Simple shear → C44 via σ6 ≈ 2 C44 ε
        def C44_from_shear(component_index, eps_list):
            # component_index: 5→xy, 4→xz, 3→yz in Voigt order [xx,yy,zz,yz,xz,xy]
            sig6 = []
            for eps in eps_list:
                Bp = np.eye(3)
                Bm = np.eye(3)
                # define symmetric shear in the chosen plane
                if component_index == 5:  # xy
                    Bp[0, 1] = Bp[1, 0] = eps
                    Bm[0, 1] = Bm[1, 0] = -eps
                elif component_index == 4:  # xz
                    Bp[0, 2] = Bp[2, 0] = eps
                    Bm[0, 2] = Bm[2, 0] = -eps
                elif component_index == 3:  # yz
                    Bp[1, 2] = Bp[2, 1] = eps
                    Bm[1, 2] = Bm[2, 1] = -eps
                sp = stress_of(C0 @ Bp)
                sm = stress_of(C0 @ Bm)
                sig6.append(((sp[component_index] - sm[component_index]) / 2.0))
            coeff = np.polyfit(eps_list, sig6, 1)[0]  # ≈ 2 C44 ε
            return 0.5 * coeff

        C44_xy = C44_from_shear(5, eps_list)
        C44_xz = C44_from_shear(4, eps_list)
        C44_yz = C44_from_shear(3, eps_list)
        C44 = float(np.mean([C44_xy, C44_xz, C44_yz]))

        # 3) Isotropic strain → bulk modulus K via P ≈ 3 K ε, P = −(σ1+σ2+σ3)/3
        P_eps = []
        for eps in eps_list:
            Bp = (1 + eps) * np.eye(3)
            Bm = (1 - eps) * np.eye(3)
            sp = stress_of(C0 @ Bp)
            sm = stress_of(C0 @ Bm)
            P = (-(np.mean(sp[:3]) - np.mean(sm[:3])) / 2.0)
            P_eps.append(P)
        K = (np.polyfit(eps_list, P_eps, 1)[0]) / 3.0  # GPa

        # 4) Recover C11 and C12 from K and (C11 - C12)
        S = 3.0 * K
        D = C11_minus_C12
        C12 = (S - D) / 3.0
        C11 = D + C12

        # 5) VRH shear modulus for cubic
        G_V = (C11 - C12 + 3.0 * C44) / 5.0
        G_R = 5.0 * (C11 - C12) * C44 / (4.0 * C44 + 3.0 * (C11 - C12))
        G = 0.5 * (G_V + G_R)

        return {
            "C11": float(C11), "C12": float(C12), "C44": float(C44),
            "K": float(K), "G": float(G),
            "C44_xy": float(C44_xy), "C44_xz": float(C44_xz), "C44_yz": float(C44_yz)
        }

    def determineCrystalStructure(self, poscar_path="POSCAR"):
        atoms = self.ideal_atoms
        [a, b, c, ang_bc, ang_ac, ang_ab] = atoms.get_cell_lengths_and_angles()
        ang_bc = round(ang_bc, 0)
        ang_ac = round(ang_ac, 0)
        ang_ab = round(ang_ab, 0)
        if (a == b and a == c and b == c) and (ang_bc == 90 and ang_ac == 90 and ang_ab == 90): # Assumes primitive cell in POSCAR
            return "bcc"
        elif (a == b and a == c and b == c) and (ang_bc == 60 and ang_ac == 60 and ang_ab == 60):   # Assumes primitive cell in POSCAR
            return "fcc"

        # Determining typ of unit cell
        """
        if (a == b and a == c and b == c) and (ang_bc == 90 and ang_ac == 90 and ang_ab == 90):
            return "cubic"
        elif (a == b and a != c and b != c) and (ang_bc == 90 and ang_ac == 90 and ang_ab == 90):
            return "tetragonal"
        elif (a != b and a != c and b != c) and (ang_bc == 90 and ang_ac == 90 and ang_ab == 90):
            return "orthorhombic"
        elif (a != b and a != c and b != c) and (ang_bc == 90 and ang_ac != 90 and ang_ab == 90):
            return "monoclinic"
        elif (a == b and a != c and b != c) and (ang_bc == 90 and ang_ac == 90 and ang_ab == 120):
            return "hexagonal"
        elif (a == b and a == c and b == c) and (ang_bc != 90 and ang_ac != 90 and ang_ab != 90):
            return "rhombohedral"
        elif (a != b and a != c and b != c) and (ang_bc != ang_ac and ang_bc != ang_ab and ang_ac != ang_ab and ang_bc != 90 and ang_ac != 90 and ang_ab != 90):
            return "triclinic"
        else:
            return "Unknown/Other"
        """

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


def ParseLogFile(log_path: str = "output.log"):
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
