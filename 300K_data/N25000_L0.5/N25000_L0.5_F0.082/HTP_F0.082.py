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

import time
import logging
import SimulationInput
import sys
import json
import numpy as np
from ASEWrappers import LangevinIntegrator, MACEPotential, AtomicStructure
from MDClasses import EquilibriumRun
from httk.external import ase_glue

def main():

    """
    Runs the high-throughput program.
    Constant parameters: number of steps, temp_k = 0, friction, time_steps.
    Saves the initial and final configurations in a cif file.
    """
    pre_time = time.time()
    # Adjustable parameters
    num_steps = 25000
    temp_k = 300
    time_steps = 0.5

    log = logging.getLogger(__name__)

    poscar_path = sys.argv[1]
    mace_path = sys.argv[2]
    try:
        lang_int = LangevinIntegrator(time_steps, temp_k, 0.082)
    except Exception as err:
        log.error(f"Integrator cannot be loaded: {err}")
        exit(1)

    mace_potential = MACEPotential(mace_path)
    settings = SimulationInput.SimulationSettings(num_steps, mace_potential, lang_int, "NVT", False)
    key, defect_index = 0, 0

    try:
        with open("key", "r", encoding="utf-8") as f:
            key = int(f.readline().rstrip("\n"))

        with open("defect_info", "r", encoding="utf-8") as f:
            defect_index = int(f.readline().rstrip("\n"))

    except Exception as err:
        log.error(f"Info files not found")

    # Load in the initial structure
    try:
        atomic_structure = AtomicStructure.fromFile(path=poscar_path, potential=mace_potential)

    except Exception as err:
        log.error(f"Cannot read the atomic structure, check if you have atomic a structure file: {err}")
        exit(1)

    # Get lattice parameter
    lattice_constant = atomic_structure.lattice_constant

    # For the atomic structure from wrapper for the initial structure
    E_pre = atomic_structure.potential_energy
    atomic_structure_atoms = atomic_structure.getAtoms()

    sorted_z_list = sorted([row[2] for row in atomic_structure_atoms.get_positions()])
    pre_factor = sorted_z_list[-1] - sorted_z_list[0]
    httk_pre = ase_glue.ase_atoms_to_structure(atomic_structure_atoms, hall_symbol="P 1")
    httk_pre.io.save("pre_structure.cif")

    # Run the simulation

    equil_MD = EquilibriumRun(settings=settings)
    equil_structure = equil_MD.run(atomic_structure, settings.num_steps, store_traj = True, check_conv=True, check_expansion=True)

    E_post = equil_structure.potential_energy

    # Calculate the expansion factor
    equil_structure_atoms = equil_structure.getAtoms()
    sorted_z_list = sorted([row[2] for row in equil_structure_atoms.get_positions()])
    post_factor = sorted_z_list[-1] - sorted_z_list[0]
    expansion_factor = post_factor / pre_factor

    # Calculate the depth of the defect
    host_array = equil_structure_atoms.get_positions()
    defect_z = host_array[defect_index][2]
    host_array = np.delete(host_array, defect_index, 0)
    sorted_z_list = sorted([row[2] for row in host_array])
    depth = (2*defect_z - (sorted_z_list[-1] + sorted_z_list[0])) / (sorted_z_list[-1] - sorted_z_list[0])

    # Save the equilibrium structure and save it in a cif file
    httk_post = ase_glue.ase_atoms_to_structure(equil_structure_atoms, hall_symbol="P 1")
    httk_post.io.save("post_structure.cif")
    with open("Output.txt",'r') as o:
        # data[0] contains the convergence criterium. 0 is energy convergence, 1 is expansion factor and 2 is time_out
        # data[1] contains the amount of steps
        data = o.read().split()
        result = {
            "MDScreenResult": {
                "key": key,
                "energy": E_post,
            },
            "MDAbadParameters": {
                "key": key,
                "depth": depth,
                "expansion_factor": expansion_factor,
                "defect_index": defect_index,
                "lattice_constant": lattice_constant,
                "convergence_criterium": int(data[0]),
                "number_of_steps": int(data[1]),
                "time": float(time.time()-pre_time),
            },
            "MDQuantities": {
                "key" : key,
                "temperature": equil_structure.temperature,
                "energy_pot": equil_structure.potential_energy,
                "energy_kinetic": equil_structure.kinetic_energy,
                "msd": equil_structure.computeMSD(atomic_structure),
                "internal_pressure": equil_structure.internal_pressure,
            }
        }
    return result
if __name__=="__main__":
    result = main()
    with open("result.json", "w", encoding="utf-8") as f:
        f.write(json.dumps(result, indent=4))
