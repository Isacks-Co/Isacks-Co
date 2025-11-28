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
from MDClasses import EquilibriumRun, SampleRun, StrecthRun
from PreProcessing import PreProcessing


def main():
    try:
        log = logging.getLogger(__name__)
        PP = PreProcessing(sys.argv)

        equil_settings, sample_settings, stretch_settings = PP.createSettings()
        atomic_structure = PP.atomic_structure

        log.info(f"Structure and settings sucessfully loaded")
    except Exception as err:
        log.error(f"Preprocessing failed: {err}")  # should probably add the err, here instead
        exit(1)
    try:
        # TODO THIS WILL GET GROUPED USING A MDMANAGER CLASS

        equil_MD = EquilibriumRun(settings=equil_settings)
        sample_MD = SampleRun(settings=sample_settings)
        stretch_MD = StrecthRun(settings=stretch_settings)

        log.info("Relaxing structure")
        equil_struct = equil_MD.run(atomic_structure, equil_settings.num_steps, init_vel=True)
        log.info("Sampling structure")
        sample_data = sample_MD.run(equil_struct, sample_settings.num_steps)
        log.info("Running stretch sequence")
        sample_data.storeTxtFile()
        C_matrix = stretch_MD.run(equil_struct)
        log.info("MD done")
        log.info(f"Stored results in {equil_struct.label}/Outputfiles")
    except Exception as err:
        log.error(f"Simulation failed: {err}")  # should probably add the err, here instead
        exit(1)

    log.info("Simulation done")


if __name__ == "__main__":
    main()
