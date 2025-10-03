import sys
from SourceCode.PreProcessing import PreProcessing
from SourceCode.MDBase import MDBase
from SourceCode.PostProcessing import PostProcessing
from SourceCode.logger import logger_setup


if __name__ == "__main__":
    log = logger_setup()


    # Clear previous simulation outputs to ensure a clean run


    if len(sys.argv[1:]) % 2 != 0:
        raise AssertionError("Missing arguments") # Maybe other exception is better

    settings = "Settings/settings.json"
    poscar = "AtomicStructure/FCC_Al.vasp"
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

    try:
        log.info("Reading settings and setting atomic structures ")
        PP = PreProcessing(settings, poscar, flags)
        PP_copy = PP
        ens_log = PP.settings.get('Ensemble', f"Equil={PP.settings.get('EquilEnsemble','?')}, Prod={PP.settings.get('ProductionEnsemble','?')}")
        log.info("Setting ensemble(s): %s and passing relevant parameters", ens_log)
        MD = PP.createMD()
    except Exception as err:
        log.error(f"Preprocessing failed: {err}") #should probably add the err, here instead
        exit(1)
    try:
        for path in (PP.settings["Output_file"] + ".traj", PP.settings["Output_file"] + ".log"):
            log.info("Creating file: %s", path)
            with open(path, "w"):
                pass
    except Exception:
        pass

    try:
        MD.runMD(PP.atoms)
    except Exception as err:
        log.error(f"Simulation failed: {err}") #should probably add the err, here instead
        exit(1)

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
        # Specific heat capacity
        PostViz.computeSpecificHeat()
    except Exception as err:
        log.error(f"computeSpecificHeat failed: {err}")
        pass


    try:
        # Debye temperature
        PostViz.computeDebyeTemperature()
    except Exception as err:
        log.error(f"computeDebyeTemperature failed: {err}")
        pass

    log.info("Simulation done")
