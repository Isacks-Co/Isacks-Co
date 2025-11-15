import sys
from pathlib import Path
import shutil
import os




    


if __name__ == "__main__":

    current_dir = Path(str(os.getcwd()))
    structure_file = Path(sys.argv[1])
    settings_file = Path(sys.argv[2])
    print(settings_file)
    
    
    current_sim = current_dir / "currentSimulation"
    
    current_sim.mkdir(exist_ok=True)

    setup_files = current_sim / "SetupFiles"
    setup_files.mkdir(exist_ok=True)
    output_files = current_sim / "OutputFiles"
    output_files.mkdir(exist_ok=True)
    # Move structure file
    shutil.copy(str(structure_file), str(setup_files / "POSCAR"))

    # Move settings file
    shutil.copy(str(settings_file), str(setup_files / "settings.json"))
