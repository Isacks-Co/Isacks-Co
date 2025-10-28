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
