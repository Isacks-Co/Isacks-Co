import json
from ase.io.vasp import read_vasp
from ase.visualize import view



class PreProcessing:
    def __init__(self, input_settings, input_structure, flags):
        self.expected_keys = {"-T": "Temperature", "-E": "Ensemble"}
        self.readSettings(input_settings)
        self.readAtomicStructure(input_structure)
        self.readTerminalInput(flags)
        self.printInput()


    def readSettings(self, input_settings):
        with open(input_settings, "r") as file:
            self.settings = json.load(file)
        for key in self.expected_keys.values():
            if not key in self.settings.keys():
                raise ValueError(f"Missing setting: {key}")
            

    def readAtomicStructure(self, input_structure):
        self.atoms = read_vasp(input_structure)
    
    
    def printInput(self):
        for key,value in self.settings.items():
            print(f"{key} : {value}")

    def readTerminalInput(self, flags):
        for i in range(0,len(flags),2):
            self.settings[self.expected_keys[flags[i]]] = flags[i+1]

        

if __name__ == "__main__":
    PreProcessing("settings.json", "poscar")
    
