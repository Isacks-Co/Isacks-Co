import sys
from PreProcessing import PreProcessing
from MDBase import MDBase
from PostProcessing import PostProcessing


if __name__ == "__main__":
    if len(sys.argv[1:]) % 2 != 0:
        raise AssertionError("Missing arguments") # Maybe other exception is better

    settings = "settings.json"
    poscar = "POSCAR"
    try:
        flags = sys.argv[1:]
    except IndexError:
        flags = None
        pass

    try:
        PP = PreProcessing(settings, poscar, flags)
        MD = PP.createMD()
    except Exception as err:
        print(f"Preprocessing failed: {err}")
        exit(1)

    try:
        MD.runMD(PP.atoms)
    except Exception as err:
        print(f"Simulation failed: {err}")
        exit(1)

    try:
        PostViz = PostProcessing('data.traj') # (TODO) Hardcoded but settings.json will contain file name
        PostViz.vizualize()
    except Exception as err:
        print(f"Postprocessing failed: {err}")
        exit(1)
