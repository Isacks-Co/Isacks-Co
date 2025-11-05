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

import sys
from PreProcessing import PreProcessing
from MDBase import MDBase
from newPostProcessing import PostProcessing
import logging


if __name__ == "__main__":
    try:
        log = logging.getLogger(__name__)
        PP = PreProcessing(sys.argv)
        settings = PP.createSettings()
        log.info(f"Settings loaded :\n{settings}")
    except Exception as err:
        log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)
    try:
        MD = MDBase(settings)
        MD.runMD(PP.atoms)
        log.info("MD done")
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)
    try:
        Postviz = PostProcessing(settings)
        
    except Exception as err:
        log.error(f"Calculating quantities failed: {err}") #should probably add the err, here instead
        exit(1)
    log.info("Simulation done")
