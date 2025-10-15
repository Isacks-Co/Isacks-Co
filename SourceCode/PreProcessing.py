import json
import sys

from ase.io import read
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
import logging
from simulationInput import NPTSettings,NVESettings,NVTSettings
from inputParser import InputParser
log = logging.getLogger(__name__)

class PreProcessing:
    """
    Class to handle all preprocessing for the MD simulation
    Idea is that this class stores all information in the two objects settings and atoms which are supposed to input
    for the MD class
    """

    def __init__(self, args):

        #Init argparser, all inputs from terminal available in self.argparser.args (dict)
        self.argparser = InputParser(args)

        #Init settings and atomic structure
        self.settings = self.readSettings(self.argparser.args["input_settings"])
        self.atoms = self.readAtomicStructure(self.argparser.args["input_structure"])

        #Physical check of the input
        self.sanityCheckAtomicStructure()
        self.sanityCheckSettings()
        self.printInput()
        

    def readSettings(self, input_settings):
        """Reads settings from json file, checks all expected settings present. Overwrite settings file if a terminal flag is set."""
        log.info("Reading settings file: %s", self.argparser.args["input_structure"])
        try:
            with open(input_settings, "r") as file:
                temp_settings = json.load(file)
        except FileNotFoundError:
            log.error("Settings file not found: %s", input_settings)
            raise FileNotFoundError(f"File {input_settings} not found, please check it exists")
        for key, value in temp_settings.items():
            if not key in self.argparser.args.keys():
                log.error("Got unexpected setting input: %s", key)
                raise ValueError(f"Got unexpected setting input: {key}")
            elif self.argparser.args[key] != None:
                temp_settings[key] = self.argparser.args[key]

        log.debug("Settings loaded: %r", temp_settings)
        return temp_settings
    

    def readAtomicStructure(self, input_structure):
        """Reads atomic structure from a file, and extend cell according to supercell setting"""

        try:
            
            log.info("Reading atomic structure from: %s", input_structure)
            atoms = read(input_structure) * tuple(self.settings["Supercells"])
            atoms.pbc = True #TODO, Hard coded pbc always true for now.
            with open(input_structure, "r") as file: # Manually read the first line and add as a comment.
                structure_name = file.readline().strip()
                atoms.info["comment"] = structure_name
            return atoms
        except FileNotFoundError:
            log.error("Structure file not found: %s", input_structure)
            raise FileNotFoundError(f"File {input_structure} not found, please check it exists")
        except Exception as err:
            error_msg = f"Atomic structure file format could not be read {err}"
            log.error(error_msg)
            raise RuntimeError(error_msg)

    def printInput(self):
        """Print out all settings to the terminal for validation"""
        for key, value in self.settings.items():
            log.info(f"{key} : {value}")
        log.info(f"Number of atoms: {len(self.atoms)}")

    def createSettings(self):
        log.info("Creating Settings object for ensemble: %s", self.settings["Ensemble"])
        match self.settings["Ensemble"]:
            case "NVE":
                
                return NVESettings(init_temp=self.settings["Temperature"], potential=self.settings["Potential"],
                                      timestep=self.settings["Timestep"], num_steps=self.settings["Number_of_steps"],
                                      interval=self.settings["Sample_interval"], output_file=self.settings["Output_file"],
                                      supercells= self.settings["Supercells"])
            case "NVT":
                
                return NVTSettings(temperature=self.settings["Temperature"], potential=self.settings["Potential"],
                                      timestep=self.settings["Timestep"], num_steps=self.settings["Number_of_steps"],
                                      interval=self.settings["Sample_interval"], output_file=self.settings["Output_file"],
                                      friction=self.settings["Friction"],
                                      supercells= self.settings["Supercells"])
            case "NPT":
                
                return NPTSettings(temperature=self.settings["Temperature"], potential=self.settings["Potential"],
                                      timestep=self.settings["Timestep"], num_steps=self.settings["Number_of_steps"],
                                      interval=self.settings["Sample_interval"], output_file=self.settings["Output_file"],
                                      pressure=self.settings["Pressure"], compressibility=self.settings["Compressibility"],
                                      supercells= self.settings["Supercells"])
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
        elif self.settings["Number_of_steps"] < 0 or not isinstance(self.settings["Number_of_steps"], int):
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

        lengths = np.array([a / i for a, i in zip(cell.lengths(), self.settings["Supercells"])])

        if np.any(angles <= 0) or np.any(angles >= 180):  # Check that lattice angles are between 0 and 180
            raise ValueError("Invalid Lattice: Lattice angles must be between 0 and 180 degrees")
        elif np.any(lengths <= 0) or np.any(lengths >= 10):  # Check so that lattice constants are not >= 10 Å (or <= 0)
            raise ValueError("Invalid Lattice: Lattice constants need to be positive and < 10")

    def checkDistances(self):
        """
        Checks that interatomic distances are reasonable. No atomic overlap
        """
        if len(self.atoms) <= 5000:  # Gets really expensive to compute interatomic distances at larger numbers
            distances_matrix = self.atoms.get_all_distances()
            upper_indeces = np.triu_indices(len(distances_matrix), k=1)
            flat_distances = distances_matrix[upper_indeces]
            if np.any(flat_distances <= 0.5): # Not sure exactly what is a reasonable threshold as atomic radius varies alot. currently 0.5 Å
                raise ValueError("Invalid atomic configuration: Atomic overlap")


if __name__ == "__main__":
    PreProcessing(sys.argv)
