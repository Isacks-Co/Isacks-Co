

from ASEWrappers.potential import EMTPotential,LennardJonesPotential
from Utils import LJParams
from ASEWrappers import VelocityVerletIntegrator,LangevinIntegrator,IsotropicMTKNPTIntegrator
from simulationInput import SimulationSettings
from Utils.inputParser import InputParser
from ASEWrappers import AtomicStructure
import json
import sys
import logging

from ase.io import read
import numpy as np


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
        print(self.settings)
        print(self.argparser.args["input_structure"])
        self.atomic_structure = self.readAtomicStructure(self.argparser.args["input_structure"])

        #Physical check of the input
        #self.sanityCheckAtomicStructure() #TODO Problems with new structure
        #self.sanityCheckSettings()


    def readSettings(self, input_settings):
        """Reads settings from json file, checks all expected settings present. Overwrite settings file if a terminal flag is set."""
        log.info("Reading settings file: %s", self.argparser.args["input_settings"])
        try:
            with open(input_settings, "r") as file:
                temp_settings = json.load(file)
        except FileNotFoundError:
            log.error("Settings file not found: %s", input_settings)
            raise FileNotFoundError(f"File {input_settings} not found, please check it exists")
        for key in temp_settings.keys():
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
            atomic_structure = AtomicStructure.fromFile(input_structure,pbc = True, supercells= self.settings["Supercells"],potential=self.getPotential())
            return atomic_structure
        
        #TODO Add this again if needed
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
        except Exception:
            error_msg = f"Atomic structure file format could not be read"
            
            log.error(error_msg)
            raise RuntimeError(error_msg)

    def getPotential(self):
        match self.settings["Potential"]:
            case "LJ":
                lj_params = LJParams(material=self.atoms.info["comment"].split()[-1]) #TODO Problem
                return LennardJonesPotential(eps = lj_params["epsilon_ev"],sigma=lj_params["sigma_A"],rc=lj_params["rc_A"],ro = lj_params["ro_A"])
            case "EMT":
                return EMTPotential()
            
        
    def getIntegrator(self,ensemble):
        match ensemble:
            case "NVE":
                return VelocityVerletIntegrator(timestep=self.settings["Timestep"])
            
            case "NVT":
                return LangevinIntegrator(timestep= self.settings["Timestep"],temperature_K=self.settings["Temperature"],friction=self.settings["Friction"])
            
            case "NPT":
                return IsotropicMTKNPTIntegrator(timestep=self.settings["Timestep"],temperature_K=self.settings["Temperature"],pressure=self.settings["Pressure"],pdamp=self.settings["Pdamp"],tdamp=self.settings["Tdamp"])
            
    def createSettings(self):
        #log.debug("Creating Settings object for ensemble: %s", self.settings["Ensemble"])
        
        potential = self.getPotential()
        self.equil_settings = SimulationSettings(num_steps=10000,potential=potential,integrator=self.getIntegrator("NPT"))
        self.sample_settings = SimulationSettings(num_steps=self.settings["Number_of_steps"],potential=potential,integrator=self.getIntegrator("NVT"))
        self.stretch_settings = SimulationSettings(num_steps=self.settings["Number_of_steps"],potential=potential,integrator=self.getIntegrator("NVT"))
        return self.equil_settings,self.sample_settings,self.stretch_settings

  



    def sanityCheckSettings(self): #TODO Put in the respective classes like integrator and atomic structure
        """
        Sanity check for the settings.json file. Makes sure that we only use EMT for
        valid metals. Also checks that relevant values are non-negative.
        """
        if self.settings["Potential"] == "EMT":
            elements = self.atoms.get_atomic_numbers()
            
            if not np.all(np.isin(elements,[13, 28, 29, 46, 47, 78, 79])): # Check if the elements are supported for EMT potential
                raise ValueError(f"Invalid potential: EMT potential only available for Al, Cu, Ag, Au, Ni, Pd, Pt.")
        if self.settings["Temperature"] > 3000:
            raise ValueError(f"Invalid temperature: Exceeds 3000K")
        elif self.settings["Temperature"] < 0:
            raise ValueError(f"Invalid temperature: Negative temperature")
        elif self.settings["Pressure"] < 0:
            raise ValueError(f"Invalid pressure: Pressure has to be non-negative")
        elif self.settings["Timestep"] < 0:
            raise ValueError(f"Invalid timestep: timestep has to be non-negative")
        elif self.settings["Friction"] < 0:
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
