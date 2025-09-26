import json
from MDBase import MDBase
from ase.io.vasp import read_vasp
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
import logging

log = logging.getLogger("MD")

class PreProcessing:
    """
    Class to handle all preprocessing for the MD simulation
    Idea is that this class stores all information in the two objects settings and atoms which are supposed to input
    for the MD class
    """

    def __init__(self, input_settings, input_structure, flags):
        self.expected_keys = {"-T": "Temperature", "-E": "Ensemble", "-P" : "Pressure",
                               "-POT" : "Potential", "-TS" : "Timestep", "-N" : "Number_of_steps",
                               "-F" : "Friction", "-C" : "Compressibility", "-I" : "Interval", "-O" : "Output_file",
                               "_A" : "Attachments", "-ES" : "EquilSteps", "-S":"Supercells"}



        log.info("Reading settings file: %s", input_settings)

        self.settings = self.readSettings(input_settings)
        print(self.settings["Supercells"])
        self.atoms = self.readAtomicStructure(input_structure) * tuple(self.settings["Supercells"])
        self.atoms.pbc = True
        self.readTerminalInput(flags)
        self.sanityCheckAtomicStructure()
        self.sanityCheckSettings()
        self.printInput()

    def readSettings(self, input_settings):
        """Reads settings from json file, checks all expected settings present"""
        try:
            with open(input_settings, "r") as file:
                temp_settings = json.load(file)
        except FileNotFoundError:
            log.error("Settings file not found: %s", input_settings)
            raise FileNotFoundError(f"File {input_settings} not found, please check it exists")
        for input in temp_settings.keys():
            if not input in self.expected_keys.values():
                log.error("Got unexpected setting input: %s", input)
                raise ValueError(f"Got unexpected setting input: {input}")

        log.debug("Settings loaded: %r", temp_settings)
        return temp_settings

    def readAtomicStructure(self, input_structure):
        """Reads atomic structure from a file with POSCAR structure"""
        try:
            log.info("Reading atomic structure from POSCAR: %s", input_structure)
            return read_vasp(input_structure)
        except FileNotFoundError:
            log.error("Structure file not found: %s", input_structure)
            raise FileNotFoundError(f"File {input_structure} not found, please check it exists")

    def printInput(self):
        """Print out all settings to the terminal for validation"""
        for key, value in self.settings.items():
            print(f"{key} : {value}")
        log.info(f"Number of atoms: {len(self.atoms)}")

    def readTerminalInput(self, flags):
        """Overwrites self.settings if other settings was received from terminal"""
        if flags:
            for i in range(0, len(flags), 2):
                if flags[i] not in self.expected_keys.keys():
                    log.error(f"Flag is invalid: {flags[i]}")
                    raise ValueError(f"Flag is invalid: {flags[i]}")
                self.settings[self.expected_keys[flags[i]]] = flags[i + 1]

    def createMD(self):
        """Init MD objects, throws errors if crucial setting is missing."""
        NVE_settings = ["Temperature", "Potential", "Timestep", "Number_of_steps", "Interval", "Output_file"]
        NVT_settings = ["Temperature", "Potential", "Timestep", "Number_of_steps", "Interval", "Output_file", "Friction"]
        NPT_settings = ["Temperature", "Potential", "Timestep", "Number_of_steps", "Interval", "Output_file", "Pressure", "Compressibility"]
        log.info("Creating MD object for ensemble: %s", self.settings["Ensemble"])
        match self.settings["Ensemble"]:
            case "NVE":
                for setting in NVE_settings:
                    if setting not in self.settings.keys():
                        raise ValueError(f"Missing the setting: {setting}")
                return MDBase.initNVE(temperature=self.settings["Temperature"], pot_str=self.settings["Potential"],
                                      timestep=self.settings["Timestep"], steps=self.settings["Number_of_steps"],
                                      interval=self.settings["Interval"], output_file=self.settings["Output_file"],
                                      equilibrium_steps = self.settings["EquilSteps"])
            case "NVT":
                for setting in NVT_settings:
                    if setting not in self.settings.keys():
                        raise ValueError(f"Missing the setting: {setting}")
                return MDBase.initNVT(temperature=self.settings["Temperature"], pot_str=self.settings["Potential"],
                                      timestep=self.settings["Timestep"], steps=self.settings["Number_of_steps"],
                                      interval=self.settings["Interval"], output_file=self.settings["Output_file"],
                                      friction=self.settings["Friction"], equilibrium_steps = self.settings["EquilSteps"])
            case "NPT":
                for setting in NPT_settings:
                    if setting not in self.settings.keys():
                        raise ValueError(f"Missing the setting: {setting}")
                return MDBase.initNPT(temperature=self.settings["Temperature"], pot_str=self.settings["Potential"],
                                      timestep=self.settings["Timestep"], steps=self.settings["Number_of_steps"],
                                      interval=self.settings["Interval"], output_file=self.settings["Output_file"],
                                      pressure_Pa=self.settings["Pressure"], compressibility=self.settings["Compressibility"],
                                      equilibrium_steps = self.settings["EquilSteps"])
            case _:
                log.error("Invalid ensemble setting: %s", self.settings["Ensemble"])
                raise ValueError(f"Invalid ensemble setting: {self.settings['Ensemble']}")

    def sanityCheckSettings(self):
        """
        Sanity check for the settings.json file. Makes sure that we only use EMT for
        valid metals. Also checks that relevant values are non-negative.
        """
        if self.settings["Potential"] == "EMT":
            elements = self.atoms.get_atomic_numbers()
            print(elements)
            if not np.all(np.isin(elements,[13, 28, 29, 46, 47, 78, 79])): # Check if the elements are supported for EMT potential
                raise ValueError(f"Invalid potential: EMT potential only availible for Al, Cu, Ag, Au, Ni, Pd, Pt.")
        if self.settings["Temperature"] > 3000:
            raise ValueError(f"Invalid temperature: Exceeds 3000K")
        elif self.settings["Temperature"] < 0:
            raise ValueError(f"Invalid temperature: Negative temperature")
        elif self.settings["Pressure"] < 0:
            raise ValueError(f"Invalid pressure: Pressure has to be non-negative")
        elif self.settings["Compressibility"] < 0:
            raise ValueError(f"Invalid compressibility: Compressibility has to be non-negative")
        elif self.settings["Friction"] < 0:
            raise ValueError(f"Invalid timestep: timestep has to be non-negative")
        elif self.settings["Timestep"] < 0:
            raise ValueError(f"Invalid friction: Friction has to be non-negative")
        elif self.settings["Number_of_steps"] < 0 or not isinstance(self.settings["Number_of_steps"],int):
            raise ValueError(f"Invalid number of steps: Has to be a positive integer")

    def sanityCheckAtomicStructure(self):
        """
        Sanity check for the input atomic structure.
        Such as valid lattice angles, constants and atomic positions
        """
        self.checkLattice()
        self.checkDistances()

    def checkLattice(self):
        """
        Check that the lattice is valid.
        """
        cell = self.atoms.get_cell()
        angles = cell.angles()

        lengths = np.array([a /i for a,i in zip(cell.lengths(), self.settings["Supercells"])])

        if np.any(angles <= 0) or np.any(angles >= 180): # Check that lattice angles are between 0 and 180
            raise ValueError("Invalid Lattice: Lattice angles must be between 0 and 180 degrees")
        elif np.any(lengths <= 0) or np.any(lengths >= 10): # Check so that lattice constants are not >= 10 Å (or <= 0)
            raise ValueError("Invalid Lattice: Lattice constants need to be positive and < 10")

    def checkDistances(self):
        """
        Checks that interatomic distances are reasonable. No atomic overlap
        """
        if len(self.atoms) <= 5000: # Gets really expensive to compute interatomic distances at larger numbers
            distances_matrix= self.atoms.get_all_distances()
            upper_indeces = np.triu_indices(len(distances_matrix), k = 1)
            flat_distances = distances_matrix[upper_indeces]
            print(flat_distances)
            if np.any(flat_distances <= 0.5): # Not sure exactly what is a reasonable threshold as atomic radius varies alot. currently 0.5 Å
                raise ValueError("Invalid atomic configuration: Atomic overlap")





if __name__ == "__main__":
    PreProcessing("settings.json", "poscar", None)
