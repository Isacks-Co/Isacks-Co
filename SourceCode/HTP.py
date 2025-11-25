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


import logging
import simulationInput
import sys
from ASEWrappers import LangevinIntegrator, MACEPotential, AtomicStructure
from MDClasses import EquilibriumRun
from httk.external import ase_glue

if __name__ == "__main__":

    """
    Runs the high-throughput program.
    Constant parameters: number of steps, temp_k = 0, friction, time_steps.
    Saves the initial and final configurations in a cif file.
    """
    # Adjustable parameters
    num_steps = 1000
    temp_k = 0
    friction = 0.1
    time_steps = 1

    log = logging.getLogger(__name__)

    poscar_path = sys.argv[1]
    mace_path = sys.argv[2]
    try: 
        lang_int = LangevinIntegrator(time_steps, temp_k, time_steps)
    except Exception as err:
        log.error(f"Integrator cannot be loaded: {err}")
        exit(1)
    mace_potential = MACEPotential(mace_path)
    settings = simulationInput.SimulationSettings(num_steps, mace_potential, lang_int)

    # Load in the initial structure
    try:
        atomic_structure = AtomicStructure.fromFile(path=poscar_path, potential=mace_potential)
    except Exception as err:
        log.error(f"Cannot read the atomic structure, check if you have atomic a structure file: {err}")
        exit(1)
    # For the atomic structure from wrapper for the initial structure
    E_pre = atomic_structure.potential_energy
    atomic_structure_atoms = atomic_structure.getAtoms()
    httk_pre = ase_glue.ase_atoms_to_structure(atomic_structure_atoms, hall_symbol="P 1")
    httk_pre.io.save("pre_structure.cif")

    # Run the simulation
    equil_MD = EquilibriumRun(settings=settings)
    equil_structure = equil_MD.run(atomic_structure, settings.num_steps)

    # Save the equilibrium structure and save it in a cif file
    E_post = equil_structure.potential_energy
    equil_structure_atoms = equil_structure.getAtoms()
    httk_post = ase_glue.ase_atoms_to_structure(equil_structure_atoms, hall_symbol="P 1")
    httk_post.io.save("post_structure.cif")

    with open("energies.txt", "w") as f:
        f.write(f"E_before  {E_pre:.10f}\n")
        f.write(f"E_after   {E_post:.10f}\n")
