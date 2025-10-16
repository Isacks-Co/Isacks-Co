def AuToGPascal(pressure_au):
    """Convert pressure from eV/Å^3 to Pa."""
    return pressure_au * 1.602176634e2

def pascalToAu(pressure_Pa):
    """Convert pressure from Pa to eV/Å^3."""
    return pressure_Pa * 6.241509074e-12

def specificHeatAuToSI(Cv_au):
    """Convert specific heat from eV/(amu·K) to J/(kg·K)."""
    return Cv_au * 9.64853321e7

def selfDiffusionCoeffAuToSI(D_au):
    return D_au * 1e-5 
