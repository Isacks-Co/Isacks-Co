import json
from MDBase import MDBase
from ase.io.vasp import read_vasp
from ase.visualize import view
from ase.lattice.cubic import FaceCenteredCubic


class PreProcessing:
    """
    Class to handle all preprocessing for the MD simulation
    Idea is that this class stores all information in the two objects settings and atoms which are supposed to input
    for the MD class
    """

    def __init__(self, input_settings, input_structure, flags):
        self.expected_keys = {"-T": "Temperature", "-E": "Ensemble"}

        self.settings = self.readSettings(input_settings)
        #self.atoms = self.readAtomicStructure(input_structure)
        #self.atoms.pbc = True
        self.atoms = FaceCenteredCubic(size=(5, 5, 5), symbol="Cu", pbc=True)

        self.readTerminalInput(flags)
        self.printInput()

    def readSettings(self, input_settings):
        """Reads settings from json file, checks all expected settings present"""
        with open(input_settings, "r") as file:
            temp_settings = json.load(file)
        for key in self.expected_keys.values():
            if not key in temp_settings.keys():
                raise ValueError(f"Missing setting: {key}")
        return temp_settings

    def readAtomicStructure(self, input_structure):
        """Reads atomic structure from a file with POSCAR structure"""
        return read_vasp(input_structure)

    def printInput(self):
        """Print out all settings to the terminal for validation"""
        for key, value in self.settings.items():
            print(f"{key} : {value}")

    def readTerminalInput(self, flags):
        """Overwrites self.settings if other settings was received from terminal"""
        for i in range(0, len(flags), 2):
            self.settings[self.expected_keys[flags[i]]] = flags[i + 1]

    def createMD(self):
        match self.settings["Ensemble"]:
            case "NVE":
                return MDBase.initNVE(temperature = self.settings["Temperature"], timestep = self.settings["Timestep_fs"],
                                       steps = self.settings["Steps"], interval = self.settings["Interval"],
                                         pot_str = self.settings["Potential"], attachments = self.settings["Attachments"] )
            case "NVT":
                return MDBase.initNVT(temperature = self.settings["Temperature"], timestep = self.settings["Timestep_fs"],
                                       steps = self.settings["Steps"], interval = self.settings["Interval"],
                                       friction = self.settings["Friction"], pot_str = self.settings["Potential"],
                                         attachments = self.settings["Attachments"] )
            case "NPT":
                return MDBase.initNPT(temperature = self.settings["Temperature"], timestep = self.settings["Timestep_fs"],
                                       steps = self.settings["Steps"], interval = self.settings["Interval"],
                                        pressure_Pa = self.settings["Pressure_Pa"],
                                      compressibility = self.settings["Compressibility"], pot_str = self.settings["Potential"],  
                                      attachments = self.settings["Attachments"])



if __name__ == "__main__":
    PreProcessing("settings.json", "poscar", None)
