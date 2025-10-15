import sys
from PreProcessing import PreProcessing
from MDBase import MDBase
from newPostProcessing import PostProcessing
import logging


if __name__ == "__main__":
    

    # Clear previous simulation outputs to ensure a clean run

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
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)

    try:
        Postviz = PostProcessing(settings)
    except Exception as err:
        log.error(f"Postprocessing failed: {err}")
        exit(1)


    log.info("Simulation done")
