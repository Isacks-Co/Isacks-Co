from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.io import read
from ase import Atoms
import re
from pathlib import Path
import numpy as np
import json


class PostProcessing:
    """
    Use MD log files to visualize and analyze results.
    """
    def __init__(self, traj_file):
        try:
            self.read_traj_file = Trajectory(traj_file)
        except FileNotFoundError:
            raise FileNotFoundError(f"Trajectory file {traj_file} not found")
        
    def vizualize(self):
        view(self.read_traj_file)

    def parse_log_file(self, log_path: str = "data.log"):
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

    def compute_cohesive_energy(self, poscar_path: str = "POSCAR", settings_path: str = "settings.json") -> float:
        """
        Compute cohesive energy per atom for the crystal in POSCAR.
        Definition: E_coh = (sum_i n_i E_i^atom − E_crystal_total) / N
        Uses the same Potential specified in settings.json (EMT or LJ).
        """
        # Read structure
        atoms = read(poscar_path)

        # Load potential type from settings
        try:
            pot_str = json.loads(Path(settings_path).read_text()).get("Potential", "EMT")
        except Exception:
            pot_str = "EMT"

        # Get calculator constructor matching MDBase behavior
        pot_str_low = str(pot_str).lower()
        calc_ctor = None
        if pot_str_low in ["emt"]:
            try:
                from asap3 import EMT as asap_EMT
                calc_ctor = asap_EMT
            except Exception:
                # Fallback to ASE EMT if ASAP3 not available
                from ase.calculators.emt import EMT as ase_EMT
                calc_ctor = ase_EMT
        elif pot_str_low in ["lj", "lennardjones", "lennard_jones", "lennardjones"]:
            from ase.calculators.lj import LennardJones
            calc_ctor = LennardJones
        else:
            # Default to EMT for unknown strings
            try:
                from asap3 import EMT as asap_EMT
                calc_ctor = asap_EMT
            except Exception:
                from ase.calculators.emt import EMT as ase_EMT
                calc_ctor = ase_EMT

        # Crystal total energy
        atoms_cryst = atoms.copy()
        atoms_cryst.calc = calc_ctor()
        E_cryst = float(atoms_cryst.get_potential_energy())

        # Isolated atom reference energies per species
        symbols = atoms_cryst.get_chemical_symbols()
        unique = {}
        for s in set(symbols):
            one = Atoms(symbols=s, positions=[[0.0, 0.0, 0.0]], cell=[20.0, 20.0, 20.0], pbc=False)
            one.calc = calc_ctor()
            unique[s] = float(one.get_potential_energy())

        # Sum over composition
        N = len(symbols)
        sum_iso = sum(unique[s] for s in symbols)

        # Cohesive energy per atom
        E_coh = (sum_iso - E_cryst) / N
        return float(E_coh)
