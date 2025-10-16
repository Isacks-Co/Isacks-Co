import sys
from SourceCode.PreProcessing import PreProcessing
from SourceCode.MDBase import MDBase

from SourceCode.quantityCalculator import QuantityCalculator
import logging

if __name__ == "__main__":
    
    try:
        #log.info("Reading settings and setting atomic structures ")
        PP = PreProcessing(sys.argv)
        settings = PP.createSettings()
        log = logging.getLogger(__name__)
        log.info("Setting ensemble: %s and passing relevant parameters", PP.settings['Ensemble'])

    except Exception as err:
        #log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)
    
    try:
        MD = MDBase(settings)
        MD.runMD(PP.atoms)
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)
    try:
        Q = QuantityCalculator(settings)
        Q.getQuantities()
    except Exception as err:
        log.error(f"Calculating quantities failed: {err}") #should probably add the err, here instead
        exit(1)
    
    log.info("Simulation done")
