import sys
from PreProcessing import PreProcessing
from MDBase import MDBase


if __name__ == "__main__":
    if len(sys.argv[1:]) % 2 != 0:
        raise AssertionError("Missing arguments") # Maybe other exception is better

    settings = "settings.json"
    poscar = "POSCAR"
    try:
        flags = sys.argv[3:]
    except IndexError:
        flags = None
        pass

    PP = PreProcessing(settings, poscar, flags)
    MD = PP.createMD()
    #print(MD.potential)
    #MD.equilibriumRun(PP.atoms)
    MD.runMD(PP.atoms)
    MD.visualizeTraj()
