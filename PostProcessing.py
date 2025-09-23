from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.io import read
from ase import Atoms
import re
from pathlib import Path
import numpy as np
import json
from ase.eos import EquationOfState
from ase.calculators.emt import EMT
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
import json
import mpmath as mp
from ase.md.analysis import DiffusionCoefficient

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

        if struct == "FCC (primitive)":
            # If these are primitive fcc vectors:
            a_from_V = (4.0 * V) ** (1.0 / 3.0)
            print("a (from primitive-cell volume) = {:.6f} Å".format(a_from_V))
            return float(a_from_V)

        elif struct == "BCC (primitive)":
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

    def ComputeMSD(self):
        traj = Trajectory(self.read_traj_file)

        N = traj[0].get_global_number_of_atoms()
        r_0 = traj[0].get_positions()
        r_n = traj[-1].get_positions()
        msd = np.sum((r_n - r_0)**2) / N
        print("Mean square displacement = ", msd, " Å^2")
        return float(msd)




    def CheckLindemannCriterion(self, poscar_path: str = "POSCAR"):
        """
        Returns True if Lindemann criterion is met, False otherwise.
        """
        traj = Trajectory(self.read_traj_file)
        atoms = traj[0]
        nl = build_neighbor_list(atoms)
        old_nn_per_atom = None
        for atoms in traj:
            nl.update(atoms)
            nn_per_atom = np.min([nearest_neighbor_distance_for_atom(i, atoms, nl) for i in range(len(atoms))])
            nn_per_atom = float(nn_per_atom)
            if old_nn_per_atom is None:
                old_nn_per_atom = nn_per_atom
            elif old_nn_per_atom is not None and (nn_per_atom >= 0) and (nn_per_atom <= old_nn_per_atom):
                old_nn_per_atom = nn_per_atom

        print("Overall nearest-neighbor distance [Å]: ", old_nn_per_atom)
        msd = self.ComputeMSD()
        L = np.sqrt(msd) / old_nn_per_atom
        print("Lindemann index: ", L)
        if L >= 0.1:
            print("Lindemann index is high.")
            return True
        else:
            print("Lindemann index is low.")
            return False

    def diffusion_coefficient(self, settings_path="settings.json", log_path="data.log", ignore_initial=0):

        traj = Trajectory(self.read_traj_file)
        atoms = traj[0]
        try:
            settings = json.loads(Path(settings_path).read_text())
            timestep_fs = float(settings["Timestep_fs"])  # fs per MD step
            interval = int(settings["Interval"])  # steps between saved frames
        except Exception as e:
            raise RuntimeError(f"Failed to read timing info from {settings_path}: {e}")

        timestep_s = timestep_fs * 1e-15  # 1 fs = 1e-15 s

        nl = build_neighbor_list(atoms)
        atom_indices, offsets = nl.get_neighbors(0)

        diff = DiffusionCoefficient(traj, timestep_s, atom_indices=atom_indices, molecule=False)
        diff.calculate(ignore_n_images=ignore_initial, number_of_segments=1)

        # diff.slopes gives slope m in: (1/(2n)) <|r(t)-r(0)|^2> = D t  => slope = D
        # Actually diff.slopes is shape (n_atom_types, n_segments, 3)
        slopes = diff.slopes  # D in each direction per type and segment
        # If only one atom type assumed, take first
        D_per_dir = slopes[0, 0, :]  # [Dx, Dy, Dz]
        D_mean = D_per_dir.mean()
        print("D_per_dir = ", D_per_dir)
        print("D_mean = ", D_mean)

        return D_per_dir, D_mean
"""
    def SelfDiffusionCoefficient(self, settings_path: str = "settings.json", fit_fraction: float = 0.5):

        # Load trajectory and timing info
        traj = Trajectory(self.read_traj_file)
        try:
            settings = json.loads(Path(settings_path).read_text())
            timestep_fs = float(settings["Timestep_fs"])  # fs per MD step
            interval = int(settings["Interval"])  # steps between saved frames
        except Exception as e:
            raise RuntimeError(f"Failed to read timing info from {settings_path}: {e}")

        dt_frame_ps = timestep_fs * interval * 1e-3  # fs → ps

        T = len(traj)
        if T < 2:
            raise ValueError("Trajectory has fewer than 2 frames; cannot compute diffusion.")

        # Initial positions and accumulator for unwrapped displacements (Å)
        cum_disp = np.zeros_like(traj[0].get_positions())

        times_ps = [0.0]
        msd_vals = [0.0]

        for k in range(1, T):
            prev = traj[k - 1]
            curr = traj[k]

            if np.any(curr.get_pbc()):
                # Minimum-image displacement in fractional coords, then map to Cartesian
                sp_prev = prev.get_scaled_positions(wrap=True)
                sp_curr = curr.get_scaled_positions(wrap=True)
                dfrac = sp_curr - sp_prev
                dfrac -= np.round(dfrac)  # nearest image
                delta = dfrac @ curr.cell.array  # (N,3) Å
            else:
                delta = curr.get_positions() - prev.get_positions()

            cum_disp += delta
            msd_k = np.mean(np.sum(cum_disp ** 2, axis=1))  # average over atoms
            msd_vals.append(float(msd_k))
            times_ps.append(k * dt_frame_ps)

        # Linear fit on the tail fraction of data
        n_pts = T
        start_idx = int((1.0 - float(fit_fraction)) * (n_pts - 1))
        start_idx = max(0, min(start_idx, n_pts - 2))
        t_fit = np.array(times_ps[start_idx:], dtype=float)
        msd_fit = np.array(msd_vals[start_idx:], dtype=float)
        if t_fit.size < 2:
            raise ValueError("Not enough points for linear fit of MSD.")

        slope, intercept = np.polyfit(t_fit, msd_fit, 1)  # MSD ≈ slope * t + intercept
        D_A2_per_ps = slope / 6.0  # Å^2/ps
        D_SI = float(D_A2_per_ps * 1e-8)  # 1 Å^2/ps = 1e-8 m^2/s

        print(f"Self-diffusion: D = {D_A2_per_ps:.6e} Å^2/ps = {D_SI:.6e} m^2/s")
        return D_SI
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
            return "FCC (primitive)"
        if all(near(x, 109.471, 2.0) for x in (al, be, ga)):
            return "BCC (primitive)"

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
