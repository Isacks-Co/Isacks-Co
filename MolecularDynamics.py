import sys
from SourceCode.PreProcessing import PreProcessing
from SourceCode.MDBase import MDBase
from SourceCode.PostProcessing import PostProcessing
from SourceCode.quantityCalculator import QuantityCalculator
import logging

if __name__ == "__main__":
    


    

    try:
        #log.info("Reading settings and setting atomic structures ")
        PP = PreProcessing(sys.argv)
        settings = PP.createSettings()
        log.info("Setting ensemble: %s and passing relevant parameters", PP.settings['Ensemble'])

    except Exception as err:
        #log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)
    log = logging.getLogger(__name__)
    try:
        MD = MDBase(settings)
        MD.runMD(PP.atoms)
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)
    Q = QuantityCalculator(settings)
    Q.computeBulkModulus()
    """

    try:
        trajectory_file = PP.settings["Output_file"] + ".traj"
        data_log_file = PP.settings["Output_file"] + ".log"
        PostViz = PostProcessing(settings, poscar, trajectory_file, data_log_file) # (TODO) Hardcoded but settings.json will contain file name
        PostViz.vizualize()
    except Exception as err:
        log.error(f"Postprocessing failed: {err}")
        exit(1)

    try:
        # Lattice constant
        PostViz.computeLatticeConstant()
    except Exception as err:
        log.error(f"computeLatticeConstant failed: {err}")
        pass

    try:
        # Cohesive energy
        PostViz.computeCohesiveEnergy()
    except Exception as err:
        log.error(f"computeCohesiveEnergy failed: {err}")
        pass

    try:
        # Bulk modulus
        PostViz.computeBulkModulus()
    except Exception as err:
        log.error(f"computeBulkModulusfailed: {err}")
        pass

    try:
        # Internal pressure
        PostViz.computeInternalPressure()
    except Exception as err:
        log.error(f"computeInternalPressure failed: {err}")
        pass

    try:
        # Mean square displacement
        PostViz.computeMSD(time=-1, reference=0)
    except Exception as err:
        log.error(f"computeMSD failed: {err}")
        pass

    try:
        # Lindemann criterion
        PostViz.computeLindemannCriterion()
    except Exception as err:
        log.error(f"computeLindemannCriterion failed: {err}")
        pass

    try:
        # Self-Diffusion Coefficient
        PostViz.computeSelfDiffusionCoefficient()
    except Exception as err:
        log.error(f"computeSelfDiffusionCoefficient failed: {err}")
        pass

    try:
        # Debye temperature
        PostViz.computeDebyeTemperature()
    except Exception as err:
        log.error(f"computeDebyeTemperature failed: {err}")
        pass
    """
    log.info("Simulation done")
