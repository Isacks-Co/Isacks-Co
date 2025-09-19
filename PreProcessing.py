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
        try:
            with open(input_settings, "r") as file:
                temp_settings = json.load(file)
        except FileNotFoundError:
            raise FileNotFoundError(f"File {input_settings} not found, please check it exists")
        for input in temp_settings.keys():
            if not input in self.expected_keys.values():
                raise ValueError(f"Got unexpected setting input: {input}")
        return temp_settings

    def readAtomicStructure(self, input_structure):
        """Reads atomic structure from a file with POSCAR structure"""
        try:
            return read_vasp(input_structure)
        except FileNotFoundError:
            raise FileNotFoundError(f"File {input_structure} not found, please check it exists")

    def printInput(self):
        """Print out all settings to the terminal for validation"""
        for key, value in self.settings.items():
            print(f"{key} : {value}")

    def readTerminalInput(self, flags):
        """Overwrites self.settings if other settings was received from terminal"""
        for i in range(0, len(flags), 2):
            if flags[i] not in self.settings.keys():
                raise ValueError(f"Flag is invalid: {flags[i]}")
            self.settings[self.expected_keys[flags[i]]] = flags[i + 1]

    def createMD(self):
        match self.settings["Ensemble"]:
            case "NVE":
                # TODO If not all settings available, raise error as such
                if False:
                    raise ValueError(f"Missing the settings: __setting__")
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

            case _:
                raise ValueError(f"Invalid ensemble setting: {self.settings['Ensemble']}")
                


if __name__ == "__main__":
    PreProcessing("settings.json", "poscar", None)
