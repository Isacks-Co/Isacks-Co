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
import sys
from PreProcessing import PreProcessing
from MDManager import MDManager


def main():
    """
    Run preprocessing and launch molecular dynamics simulations.

    This function initializes preprocessing based on command-line arguments,
    prepares simulation settings, and executes the simulations through the
    MDManager interface.

    Notes
    -----
    - Command-line arguments are passed directly to the PreProcessing class.
    - Any failure during preprocessing or simulation execution results in a
      logged exception and program termination.
    """

    log = logging.getLogger(__name__)

    # --------------------
    # Preprocessing stage
    # --------------------
    try:
        # Parse input arguments and initialize preprocessing
        PP = PreProcessing(sys.argv)

        # Log quantities requested for computation
        log.info(
            f"Quantities to compute:\n {PP.settings['Compute_quantities']}\n"
        )

        # Generate simulation settings list
        settings_list = PP.createSimulationList()

        # Retrieve atomic structure prepared during preprocessing
        atomic_structure = PP.atomic_structure

        # Initialize MD manager with prepared settings
        SimulationSetup = MDManager(settings_list, atomic_structure)

    except Exception as err:
        log.exception(f"Preprocessing failed: {err}")
        exit(1)

    # --------------------
    # Simulation execution
    # --------------------
    try:
        # Run simulations and store trajectories
        SimulationSetup.run(store_traj=True)

    except Exception as err:
        log.exception(f"Simulation failed: {err}")
        exit(1)

    log.info("Simulation done")


if __name__ == "__main__":
    main()