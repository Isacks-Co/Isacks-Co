# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import json
import logging
import numpy as np
import sys
from ASEWrappers import AtomicStructure
from ASEWrappers import VelocityVerletIntegrator, LangevinIntegrator, BerendsenNPTIntegrator
from ASEWrappers.potential import EMTPotential, LennardJonesPotential, MACEPotential
from Utils import LJParams
from Utils.inputParser import InputParser
from ase.io import read
from SimulationInput import SimulationSettings

log = logging.getLogger(__name__)


class PreProcessing:
    """
    Class to handle all preprocessing for the MD simulation
    Idea is that this class stores all information in the two objects settings and atoms which are supposed to input
    for the MD class
    """

    def __init__(self, args):

        # Init argparser, all inputs from terminal available in self.argparser.args (dict)
        self.argparser = InputParser(args)

        # Init settings and atomic structure

        self.settings = self.readSettings(self.argparser.args["input_settings"])

        self.atomic_structure = self.readAtomicStructure("../SetupFiles/POSCAR")

        # Physical check of the input
        self.sanityCheckAtomicStructure()
        self.sanityCheckSettings()

    def readSettings(self, input_settings):
        """Reads settings from json file, checks all expected settings present. Overwrite settings file if a terminal flag is set."""
        log.debug("Reading settings file: %s", self.argparser.args["input_settings"])
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
            atomic_structure = AtomicStructure.fromFile(input_structure, pbc=self.settings["Simulations_config"]["PBC"],
                                                        supercells=self.settings["Simulations_config"]["Supercells"],
                                                        potential=self.getPotential())
            return atomic_structure

        except FileNotFoundError:
            log.error("Structure file not found: %s", input_structure)
            raise FileNotFoundError(f"File {input_structure} not found, please check it exists")
        except Exception:
            error_msg = f"Atomic structure file format could not be read"

            log.error(error_msg)
            raise RuntimeError(error_msg)

    def getPotential(self):  # TODO FIX
        match self.settings["Simulations_config"]["Potential"]["Kind"]:
            case "LJ":
                atoms = read("../SetupFiles/POSCAR")

                atomic_num = [(atoms.get_atomic_numbers()[0])]
                atomic_symbols = atoms.get_chemical_symbols()
                epsilon = None
                sigma = None
                rc = None
                ro = None
                if "Material" in self.settings["Simulations_config"]["Potential"]["Parameters"]:
                    material = self.settings["Simulations_config"]["Potential"]["Parameters"]["Material"]
                else:
                    material = sorted(set(atomic_symbols))[0]

                if "epsilon_eV" in self.settings["Simulations_config"]["Potential"]["Parameters"]:
                    epsilon = self.settings["Simulations_config"]["Potential"]["Parameters"]["epsilon_eV"]
                if "sigma" in self.settings["Simulations_config"]["Potential"]["Parameters"]:
                    sigma = self.settings["Simulations_config"]["Potential"]["Parameters"]["sigma"]
                if "RC" in self.settings["Simulations_config"]["Potential"]["Parameters"]:
                    rc = self.settings["Simulations_config"]["Potential"]["Parameters"]["RC"]
                if "RO" in self.settings["Simulations_config"]["Potential"]["Parameters"]:
                    ro = self.settings["Simulations_config"]["Potential"]["Parameters"]["RO"]

                lj_params = LJParams(material=material, epsilon_eV=epsilon, sigma_A=sigma, rc_A=rc, ro_A=ro)
                return LennardJonesPotential(atomic_numbers=atomic_num, epsilons=[lj_params["epsilon_eV"]],
                                             sigmas=[lj_params["sigma_A"]], rc=lj_params["rc_A"])
            case "EMT":
                return EMTPotential()

            case "MACE":
                return MACEPotential(model_path=self.settings["Simulations_config"]["Potential"]["Parameters"]["Path"])

    def getIntegrator(self, ensemble):
        match ensemble:
            case "NVE":
                return VelocityVerletIntegrator(timestep=self.settings["Simulations_config"]["Timestep"])

            case "NVT":
                return LangevinIntegrator(timestep=self.settings["Simulations_config"]["Timestep"],
                                          temperature_K=self.settings["Physical_environment"]["Temperature"],
                                          friction=self.settings["Simulations_config"]["Friction"])

            case "NPT":
                return BerendsenNPTIntegrator(timestep = self.settings["Simulations_config"]["Timestep"],
                                              temperature_K = self.settings["Physical_environment"]["Temperature"],
                                              pressure = self.settings["Physical_environment"]["Pressure"],
                                              compressibility = self.settings["Simulations_config"]["Compressibility"]
                )

    def createSimulationList(self):
        """
        Function that creates a simulation list wherer each element dictates a specific run and settings for that run.

        Returns:
            simulation_list (list) : A list of String describing simulation and SimulationSettings objects
        """
        simulation_list = []
        potential = self.getPotential()
        NEED_STRETCH = ["Moduli", "Debye"]

        self.npt_settings = SimulationSettings(num_steps=10000, potential=potential,
                                                 integrator=self.getIntegrator("NPT"), ensemble = "NPT", sample_nn = False)

        nn_sample_necessary = True if "L_crit" in self.settings["Compute_quantities"] else False
        self.nvt_settings = SimulationSettings(num_steps=self.settings["Simulations_config"]["Number_of_steps"], potential=potential,
                                                  integrator=self.getIntegrator("NVT"), ensemble = "NVT", sample_nn = nn_sample_necessary)

        if self.settings["Find_equilibrium"]:
            simulation_list.append(["Equilibrium", self.npt_settings])


        stretch_flag, sample_flag = False, False
        for element in self.settings["Compute_quantities"]:
            if element in NEED_STRETCH:
                stretch_flag = True
            else:
                sample_flag = True

        if sample_flag:
            simulation_list.append(["Sample", self.nvt_settings])

        if stretch_flag:
            simulation_list.append(["Stretch" , self.nvt_settings])

        log.info("Simulation steps that will run:")
        for simulation in simulation_list:
            log.info(f"        ╰┈➤     {simulation[0]} \n")

        return simulation_list

    def sanityCheckSettings(self):
        """
        Sanity check for the settings.json file. Makes sure that we only use EMT for
        valid metals. Also checks that relevant values is given.
        """
        if self.settings["Simulations_config"]["Potential"]["Kind"] == "EMT":
            elements = self.atomic_structure.get_atomic_numbers

            if not np.all(np.isin(elements, [13, 28, 29, 46, 47, 78,
                                             79])):  # Check if the elements are supported for EMT potential
                raise ValueError(f"Invalid potential: EMT potential only available for Al, Cu, Ag, Au, Ni, Pd, Pt.")
        if self.settings["Physical_environment"]["Temperature"] > 3000:
            raise ValueError(f"Invalid temperature: Exceeds 3000K")
        elif self.settings["Physical_environment"]["Temperature"] < 0:
            raise ValueError(f"Invalid temperature: Negative temperature")
        elif self.settings["Physical_environment"]["Pressure"] < 0:
            raise ValueError(f"Invalid pressure: Pressure has to be non-negative")
        elif self.settings["Simulations_config"]["Timestep"] < 0:
            raise ValueError(f"Invalid timestep: timestep has to be non-negative")
        elif self.settings["Simulations_config"]["Friction"] < 0:
            raise ValueError(f"Invalid friction: Friction has to be non-negative")
        elif self.settings["Simulations_config"]["Number_of_steps"] < 0 or not isinstance(self.settings["Simulations_config"]["Number_of_steps"], int):
            raise ValueError(f"Invalid number of steps: Has to be a positive integer")
        elif self.settings["Simulations_config"]["Number_of_steps"] <= 100:
            raise ValueError(f"Invalid number of steps: Has to sample over more than 100 steps")

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
        cell = self.atomic_structure.get_cell
        angles = cell.angles()

        lengths = np.array([a / i for a, i in zip(cell.lengths(), self.settings["Simulations_config"]["Supercells"])])

        if np.any(angles <= 0) or np.any(angles >= 180):  # Check that lattice angles are between 0 and 180
            raise ValueError("Invalid Lattice: Lattice angles must be between 0 and 180 degrees")
        elif np.any(lengths <= 0) or np.any(lengths >= 10):  # Check so that lattice constants are not >= 10 Å (or <= 0)
            raise ValueError("Invalid Lattice: Lattice constants need to be positive and < 10")

    def checkDistances(self):
        """
        Checks that interatomic distances are reasonable. No atomic overlap
        """
        if len(self.atomic_structure) <= 5000:  # Gets really expensive to compute interatomic distances at larger numbers
            distances_matrix = self.atomic_structure.get_all_distances
            upper_indeces = np.triu_indices(len(distances_matrix), k=1)
            flat_distances = distances_matrix[upper_indeces]
            if np.any(
                    flat_distances <= 0.5):  # Not sure exactly what is a reasonable threshold as atomic radius varies alot. currently 0.5 Å
                raise ValueError("Invalid atomic configuration: Atomic overlap")


if __name__ == "__main__":
    PreProcessing(sys.argv)
