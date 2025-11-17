import sys

from PreProcessing import PreProcessing
<<<<<<< HEAD
from MDBase import MDBase
<<<<<<< HEAD
from PostProcessing import PostProcessing
=======
from SourceCode.PostProcessing import PostProcessing
>>>>>>> ae14402 ( Moved some of the writing of quantuities to psotprocessing and changed QC to only contain functions)
=======
from MDClasses import EquilibriumRun,SampleRun,StrecthRun



import logging

def main():
    try:
        log = logging.getLogger(__name__)
        PP = PreProcessing(sys.argv)
      
        equil_settings,sample_settings,stretch_settings = PP.createSettings()
        atomic_structure = PP.atomic_structure
        
        log.info(f"Structure and settings sucessfully loaded")
    except Exception as err:
        log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)
    try:
        #TODO THIS WILL GET GROUPED USING A MDMANAGER CLASS
       

        
        
        
        
        equil_MD = EquilibriumRun(settings = equil_settings)
        sample_MD = SampleRun(settings = sample_settings)
        stretch_MD = StrecthRun(settings = stretch_settings)
        
        log.info("Relaxing structure")
        equil_struct = equil_MD.run(atomic_structure,equil_settings.num_steps) 
        log.info("Sampling structure")
        sample_data = sample_MD.run(equil_struct,sample_settings.num_steps)
        log.info("Running stretch sequence")
        sample_data.storeTxtFile()
        C_matrix = stretch_MD.run(equil_struct)
        log.info("MD done")
        log.info(f"Stored results in {equil_struct.label}/Outputfiles")
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)
    

    log.info("Simulation done")

if __name__ == "__main__":
    main()
