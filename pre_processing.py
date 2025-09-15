import json
from ase.io.vasp import read_vasp
from ase.visualize import view

class PreProcessing:
    def __init__(self, input_settings, input_structure):
        self.readSettings(input_settings)
        self.readAtomicStructure(input_structure)
        self.printInput()


    def readSettings(self, input_settings):
        with open(input_settings, "r") as file:
            self.settings = json.load(file)

    def readAtomicStructure(self, input_structure):
        self.atoms = read_vasp(input_structure)
        view(self.atoms)
    
    
    def printInput(self):
        for key,value in self.settings.items():
            print(f"{key} : {value}")

    def readTerminalInput(self):
        pass

if __name__ == "__main__":
    PreProcessing("settings.json", "poscar")
    
