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

EV_PER_A3_TO_GPA = 160.21766208

class PostProcessing:
    """
        Uses the log files generates from the Molecular Dynamics class to vizualize and analyze the data
        input trajectory file as traj_file 
    """
    def __init__(self, traj_file):
        try:
            self.read_traj_file = Trajectory(traj_file)
        except FileNotFoundError:
            raise FileNotFoundError(f"Trajectory file {traj_file} not found")
        
    def vizualize(self):
        view(self.read_traj_file)

    def ComputeCohesiveEnergy(self, poscar_path: str = "POSCAR", settings_path: str = "settings.json") -> float:
        """
        Compute cohesive energy per atom for the crystal in POSCAR.
        Definition: E_coh = (sum_i n_i E_i^atom − E_crystal_total) / N
        Note: Uses a fixed EMT calculator; does not read settings or apply fallbacks.
        """
        # Read structure
        atoms = read(poscar_path)

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

    def ComputeLatticeConstant(self, poscar_path: str = "POSCAR") -> float:
        """
        Calculates the lattice constant of a FCC crystal.
        """
        atoms = read(poscar_path)
        # 3x3 lattice matrix (rows are the lattice vectors in Cartesian Å)
        cell = atoms.cell.array  # numpy.ndarray of shape (3, 3)

        a1, a2, a3 = cell[0], cell[1], cell[2]
        # compute volume
        V = abs(np.dot(a1, np.cross(a2, a3)))

        struct = DetermineCrystalStructure(poscar_path)
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

    def ComputeBulkModulus(self, poscar_path: str = "POSCAR", scales = np.linspace(0.96, 1.04, 7)):
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
        print("Bulk modulus = {:.6f} GPa".format(bulk_modulus))
        return float(bulk_modulus)

    def ComputeInternalPressure(self):
        """
        For NVT ensemble, computes internal pressure.
        Instantaneous: P(t)  = (1/3V)[2NkT(t) + SUM_i(r_i*f_i)]
        Average:       P = (1/M)SUM_i(P(nΔt))
        """
        traj = Trajectory(self.read_traj_file)
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

    def ComputeMSD(self, time = -1, reference = None, flags = [], log_path = "data.log", settings_path = "settings.json"):
        traj = Trajectory(self.read_traj_file)
        data_log = ParseLogFile(log_path)
        if reference is None:
            ideal_atoms = PreProcessing.PreProcessing(input_settings=settings_path, input_structure=None, flags=flags)
            N = (ideal_atoms.atoms).get_global_number_of_atoms()
            r_0 = ideal_atoms.atoms.get_positions()

        else:
            N = traj[reference].get_global_number_of_atoms()
            r_0 = traj[reference].get_positions()

        r_n = traj[time].get_positions()
        msd = np.sum((r_n - r_0)**2) / N        #Å^2
        print("MSD:---------------------- ", msd)
        return msd

        """
        # Convert ASE trajectory to a multi-frame XYZ
        images = read("data.traj", index=":")
        write("data.xyz", traj)     # multi-frame XYZ
        u = mda.Universe("data.xyz")
        print("hej")
        EMSD = mda_msd.EinsteinMSD(u, select="all", msd_type="xyz", fft=True)
        EMSD.run()
        msd = EMSD.results.timeseries
        print("MSD values:---------------------- ", msd[time])
        return msd[time]
        """

    def CheckLindemannCriterion(self, poscar_path: str = "POSCAR", settings_path = "settings.json", flags = []):
        """
        Returns True if Lindemann criterion is met, False otherwise.
        """
        #atoms_initial = read(poscar_path)
        #supercell = atoms_initial.repeat([2,2,2])
        ideal_atoms = PreProcessing.PreProcessing(input_settings=settings_path, input_structure=None, flags=flags)
        supercell = ideal_atoms.atoms
        cutoffs = natural_cutoffs(supercell)
        nl_initial = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl_initial.update(supercell)
        min_dist = np.min([nearest_neighbor_distance_for_atom(j, supercell, nl_initial) for j in range(supercell.get_global_number_of_atoms())])
        min_dist = float(min_dist)
        print("First min dist-------------- ", min_dist)

        traj = Trajectory(self.read_traj_file)
        nl = build_neighbor_list(traj[0])
        for atoms in traj:
            nl.update(atoms)
            nn_per_atom = np.min([nearest_neighbor_distance_for_atom(i, atoms, nl) for i in range(len(atoms))])
            nn_per_atom = float(nn_per_atom)
            if min_dist is None:
                min_dist = nn_per_atom
            elif min_dist is not None and (nn_per_atom >= 0) and (nn_per_atom <= min_dist):
                min_dist = nn_per_atom

        print("Overall nearest-neighbor distance [Å]: ", min_dist)
        msd = self.ComputeMSD()
        l = np.sqrt(msd) / min_dist
        print("Lindemann index: ", l)
        if l >= 0.1:
            print("Lindemann index is high.")
            return True
        else:
            print("Lindemann index is low.")
            return False

    def SelfDiffusionCoefficient(self, settings_path="settings.json", log_path="data.log"):

        traj = Trajectory(self.read_traj_file)
        data_log = ParseLogFile(log_path)
        try:
            settings = json.loads(Path(settings_path).read_text())
            timestep_fs = float(settings["Timestep_fs"])  # fs per MD step
            interval = int(settings["Interval"])  # steps between saved frames
        except Exception as e:
            raise RuntimeError(f"Failed to read timing info from {settings_path}: {e}")

        timestep_list = []
        for i in range(len(data_log["Time[ps]"])):
            if round(data_log["T[K]"][i], 0) == settings["Temperature"]:
                timestep = data_log["Time[ps]"][i]
                timestep_list.append([timestep, i])

        print("Timestep list: ", timestep_list)
        print(timestep_list[0][0], "--------", timestep_list[-1][0])


        #halfway = int(round((len(timestep_list) - 1) / 2, 0))
        msd_list = [self.ComputeMSD(time=timestep_list[0][1]),
                    self.ComputeMSD(time=timestep_list[-1][1]),]
        slope = (msd_list[-1] - msd_list[0])/ (timestep_list[-1][0] - timestep_list[0][0])
        slope_m2_per_s = slope * 1e-8 / 6
        print("Self-Diffusion Coefficient------------ [Å^2/s^-1]: ", slope_m2_per_s)
        return slope_m2_per_s


        """
        images = read("data.traj", index=":")
        write("data.xyz", traj)     # multi-frame XYZ
        u = mda.Universe("data.xyz")

        u = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)
        MSD = mda_msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
        MSD.run()
        msd = MSD.results.timeseries

        nframes = MSD.n_frames
        timestep = timestep_fs * interval  # this needs to be the actual time between frames
        lagtimes = np.arange(nframes) * timestep  # make the lag-time axis
        fig = plt.figure()
        ax = plt.axes()
        # plot the actual MSD
        ax.plot(lagtimes, msd, color="black", ls="-", label=r'3D random walk')
        exact = lagtimes * 6
        # plot the exact result
        ax.plot(lagtimes, exact, color="black", ls="--", label=r'$y=2 D\tau$')
        #plt.show()
        plt.loglog(lagtimes, msd)
        #plt.show()

        #start_time = 20
        start_index = 1 #int(start_time / timestep)
        #end_time = 60
        end_index = -1 #int(end_time / timestep)
        linear_model = linregress(lagtimes[start_index:end_index],
                                  msd[start_index:end_index])
        slope = linear_model.slope
        error = linear_model.stderr
        # dim_fac is 3 as we computed a 3D msd with 'xyz'
        D = slope * 1 * 1e-8 / (2 * MSD.dim_fac)
        print(D)
        """


def ParseLogFile(log_path: str = "data.log"):
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

def DetermineCrystalStructure(poscar_path="POSCAR"):
    atoms = read(poscar_path)
    cell = atoms.cell.array
    fracs = atoms.get_scaled_positions(wrap=True)

    def lengths_angles(M):
        a, b, c = [np.linalg.norm(v) for v in M]
        def ang(u, v):
            cos = np.clip(np.dot(u, v) / (np.linalg.norm(u)*np.linalg.norm(v)), -1, 1)
            return np.degrees(np.arccos(cos))
        alpha, beta, gamma = ang(M[1], M[2]), ang(M[0], M[2]), ang(M[0], M[1])
        return (a, b, c), (alpha, beta, gamma)

    (a, b, c), (al, be, ga) = lengths_angles(cell)
    eq_len = max(abs(a-b), abs(b-c), abs(a-c)) / max(a, b, c) < 1e-3    #- Angle tolerances are ±2 degrees, fractional-position tolerance is ~1e-3 in fractional coordinates.
    near = lambda x, y, tol: abs(x - y) < tol

    # Checks if primitive fcc or primitive bcc
    if eq_len:
        if all(near(x, 60.0, 2.0) for x in (al, be, ga)):
            return "fcc"
        if all(near(x, 109.471, 2.0) for x in (al, be, ga)):
            return "bcc"

    # Conventional cubic check
    if eq_len and all(near(x, 90.0, 2.0) for x in (al, be, ga)):
        key = lambda p: tuple(np.round(np.mod(p, 1.0)/1e-3).astype(int))
        pos = {key(p) for p in fracs}
        sc  = {key([0,0,0])}
        bcc = {key([0,0,0]), key([0.5,0.5,0.5])}
        fcc = {key([0,0,0]), key([0,0.5,0.5]), key([0.5,0,0.5]), key([0.5,0.5,0])}
        if pos == fcc:
            return "FCC (conventional)"
        if pos == bcc:
            return "BCC (conventional)"
        if pos == sc:
            return "SC (conventional)"

    return "Unknown/Other"

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


def numeric_limit_infinity(f, direction='+', t0=1.0, growth=2.0, steps=8):
    """Estimate lim_{x→±∞} f(x) by sampling x_k growing geometrically.
    direction: '+' for +∞, '-' for -∞
    t0 > 0 sets the starting magnitude; growth > 1 controls how fast x grows.
    """
    import mpmath as mp

    if growth <= 1:
        raise ValueError("growth must be > 1")

    xs = []
    x = float(t0)
    for _ in range(steps):
        xs.append(x if direction == '+' else -x)
        x *= growth

    vals = [f(xx) for xx in xs]
    vals = [v for v in vals if mp.isfinite(v)]
    return (mp.mpf(sum(vals)) / len(vals)) if vals else mp.nan
