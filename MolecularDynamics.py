import sys
from PreProcessing import PreProcessing
from MDBase import MDBase
from PostProcessing import PostProcessing
import logging
import logger


if __name__ == "__main__":
    log = logger.logger_setup()


    # Clear previous simulation outputs to ensure a clean run
    try:
        for path in ("data.log", "data.traj"):
            with open(path, "w"):
                pass
    except Exception:
        pass


    if len(sys.argv[1:]) % 2 != 0:
        raise AssertionError("Missing arguments") # Maybe other exception is better

    settings = "settings.json"
    poscar = "POSCAR"
    try:
        flags = sys.argv[1:]
        if flags[0] == "--help":
            print("To run this program use : python3 MolecularDynamics.py\n")
            print("Make sure to have a settings.json and a POSCAR file in the same directory as the MD program.\n")
            print("Instead of a settings file, the program can be run as a terminal program. Values can then be assigned as terminal input. See flags below.\n")
            available_flags = {"-T": "Temperature as a float (K)", "-E": "Ensemble as a string (NVE, NPT, NPT)",
                               "-P" : "Pressure as a float (Pa)",
                               "-POT" : "Potential as a string (EMT, LJ, MACE)",
                               "-TS" : "Timestep as a float (fs)", "-N" : "Number of steps as an integer",
                               "-F" : "Friction as a float (unit?)",
                               "-C" : "Compressibility as a float (unit?)",
                               "-I" : "Interval as an integer in which information is sampled (10 means that for e.g energy is sampled every 10 steps)",
                               "-O" : "Output file, path to desired destination for data as a string"}
            print(f"Available flags :\n") # Not exactly sure what to print here more than flags.
            for key, value in available_flags.items():
                print(f"{key} : {value}")
            sys.exit(1)
        elif len(flags) % 2 != 0:
            raise AssertionError("Invalid input") # Maybe other exception is better
    except AssertionError:
        print("Invalid input, use --help for more information")
        sys.exit(1)
    except IndexError:
        flags = None
        pass

    settings = "settings.json"
    poscar = "POSCAR"

    try:
        log.info("Reading settings and setting atomic structures ")
        PP = PreProcessing(settings, poscar, flags)
        log.info("Setting ensemble: %s and passing relevant parameters", PP.settings['Ensemble'])
        MD = PP.createMD()
    except Exception as err:
        log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)

    try:
        MD.runMD(PP.atoms)
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)

    try:
        output_str = PP.settings["Output_file"] + ".traj"
        PostViz = PostProcessing(output_str) # (TODO) Hardcoded but settings.json will contain file name
        PostViz.vizualize()
        """
        # Lattice constant
        PostViz.ComputeLatticeConstant()
        # Cohesive energy
        PostViz.ComputeCohesiveEnergy()
        # Bulk modulus
        PostViz.ComputeBulkModulus()
        # Internal pressure
        PostViz.ComputeInternalPressure()
        #Mean square displacement
        PostViz.ComputeMSD(reference=0, flags=flags)
        #Lindemann index
        PostViz.CheckLindemannCriterion(flags=flags)
        #Self-Diffusion Coefficient
        PostViz.SelfDiffusionCoefficient()
        """
        #Debye temperature
        PostViz.ComputeDebyeTemperature()
    except Exception as err:
        log.error(f"Postprocessing failed: {err}")
        exit(1)

    log.info("Simulation done")
