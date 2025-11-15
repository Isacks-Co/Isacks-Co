import sys
from pathlib import Path
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
from PostProcessing import PostProcessing


import logging

def main():
    try:
        log = logging.getLogger(__name__)
        PP = PreProcessing(sys.argv)
        settings = PP.createSettings()
        atomic_structure = PP.atomic_structure
        
        log.info(f"Settings loaded :\n{settings}")
    except Exception as err:
        log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)
    try:
        #TODO THIS WILL GET GROUPED USING A MDMANAGER CLASS
       

        
        

        equil_MD = EquilibriumRun(settings = settings)
        sample_MD = SampleRun(settings = settings)
        stretch_MD = StrecthRun(settings = settings)
        
        equil_struct = equil_MD.run(atomic_structure,settings.num_steps) # Equilibrium run
        print(equil_struct.label)
        sample_data = sample_MD.run(equil_struct,settings.num_steps)
        sample_data.storeTxtFile()
        C_matrix = stretch_MD.run(equil_struct,settings.num_steps)
        log.info("MD done")
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)
    

    log.info("Simulation done")

if __name__ == "__main__":
    main()
