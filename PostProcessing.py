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
        atoms = read("POSCAR")
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
            return a_from_V

        elif struct == "BCC (primitive)":
            # If these are primitive bcc vectors:
            a_from_V = (2.0 * V) ** (1.0 / 3.0)
            print("a (from primitive-cell volume) = {:.6f} Å".format(a_from_V))
            return a_from_V

        else:
            # If these were conventional cubic vectors (orthogonal, equal length):
            a_conv = np.linalg.norm(a1)  # would equal the conventional a
            print("a (if conventional cell) = {:.6f} Å".format(a_conv))
            #Added the section immediately above, because I don't remember much of solid state physics, and want to be sure I'm doing this right.
            return a_conv

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
        return bulk_modulus

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
            number_of_atoms = len(atoms)
            volume = atoms.get_volume()
            internal_pressure.append(((1/(3*volume))*((2*number_of_atoms*e_kin) + sum(atoms.get_forces()*atoms.get_positions())))*EV_PER_A3_TO_GPA)
        internal_pressure = np.array(internal_pressure)
        average_internal_pressure = np.average(internal_pressure)
        print("Average pressure = ", average_internal_pressure, "GPa")
        return average_internal_pressure


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
