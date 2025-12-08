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


import os
import shutil
import sys
import json
from pathlib import Path

if __name__ == "__main__":
    current_dir = Path(str(os.getcwd()))
    settings_file = Path(sys.argv[1])
    with open(settings_file, "r") as file:
        settings_dict = json.load(file)
    structure_file = Path(settings_dict["Input_structure"])
    current_sim = current_dir / "currentSimulation"

    current_sim.mkdir(exist_ok=True)

    setup_files = current_sim / "SetupFiles"
    setup_files.mkdir(exist_ok=True)
    output_files = current_sim / "OutputFiles"
    output_files.mkdir(exist_ok=True)
    # Move structure file
    print(structure_file)
    print(current_dir)
    print(f"{str(setup_files)}/atomic_structure.{str(structure_file).rsplit('.', 1)[1]}")
    shutil.copy(f"{str(structure_file)}", f"{str(setup_files)}/atomic_structure.{str(structure_file).rsplit('.', 1)[1]}")



    # Move settings file
    shutil.copy(str(settings_file), str(setup_files / "settings.json"))
