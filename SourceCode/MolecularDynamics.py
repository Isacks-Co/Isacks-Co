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
    log = logging.getLogger(__name__)
    try:
        PP = PreProcessing(sys.argv)

        settings_list = PP.createSettings()
        log.info(f"Settings created")
        quantitites_to_compute = PP.settings["Compute_quantities"]
        log.info(f"Quantities to compute: {quantitites_to_compute}")
        SimulationSetup = MDManager(settings_list, quantitites_to_compute)
        log.info(f"Simulation is ready to run")
        log.info(f"Order of operations for this run: {SimulationSetup.order_of_operations}")
        log.info(f"Quantities by variation: {SimulationSetup.quants_by_variation}")

        atomic_structure = PP.atomic_structure
    except Exception as err:
        log.error(f"Preprocessing failed: {err}")  # should probably add the err, here instead
        exit(1)
    try:
        SimulationSetup.run(atomic_structure, init_vel=False, store_traj=False)
    except Exception as err:
        log.error(f"Simulation failed: {err}")  # should probably add the err, here instead
        exit(1)

    log.info("Simulation done")


if __name__ == "__main__":
    main()
