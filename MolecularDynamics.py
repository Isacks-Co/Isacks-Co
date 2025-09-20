import sys
from PreProcessing import PreProcessing
from MDBase import MDBase
from PostProcessing import PostProcessing
import logging


if __name__ == "__main__":
    logging.basicConfig( level=logging.INFO,
        format="%(levelname)s: %(message)s")

    log = logging.getLogger("MD")


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
        log.info(f"Setting the {PP.settings["Ensemble"]} and passing relevant parameters ")
        MD = PP.createMD()
    except Exception as err:
        print(f"Preprocessing failed: {err}")
        log.exception("Preprocessing failed") #should probably add the err, here instead
        exit(1)

    try:
        log.info("Starting the Molecular dynamic simulation ")
        MD.runMD(PP.atoms)
    except Exception as err:
        print(f"Simulation failed: {err}")
        log.exception("Simulation failed") #should probably add the err, here instead
        exit(1)

    try:
        PostViz = PostProcessing('data.traj') # (TODO) Hardcoded but settings.json will contain file name
        PostViz.vizualize()
    except Exception as err:
        print(f"Postprocessing failed: {err}")
        exit(1)

log.info("Simulation done")
