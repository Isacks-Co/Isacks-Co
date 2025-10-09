LJ_DB = {
    "ar": {"epsilon_eV": 0.0103, "sigma": 3.405},
    "kr": {"epsilon_eV": 0.01475, "sigma": 3.650},
    "xe": {"epsilon_eV": 0.01903, "sigma": 4.100},
    "ne": {"epsilon_eV": 0.00307, "sigma": 2.789},
    "cu": {"epsilon_eV": 0.40,  "sigma": 2.28},
    "fe": {"epsilon_eV": 0.2007, "sigma": 2.4193},
    "ni": {"epsilon_eV": 0.1729, "sigma": 1.5808},
}

def LJParams(material: str = None, *, epsilon_eV=None, sigma_A=None, rc_A=None, ro_A=None):

    if material:
        base = LJ_DB.get(material.lower())
        if not base:
            raise ValueError(f"No LJ-parametrar registered for '{material}'."
                             f"Change to any of these: "
                             f"{", ".join(LJ_DB.keys())}")
    else:
        base = {}

    eps = epsilon_eV if epsilon_eV is not None else base.get("epsilon_eV")
    sig = sigma_A   if sigma_A   is not None else base.get("sigma")
    if eps is None or sig is None:
        raise ValueError("Lennard-Jones needs epsilon_eV and sigma.")

    rc = rc_A if rc_A is not None else 2.5*float(sig)
    ro = ro_A if ro_A is not None else 0.9*float(rc)

    if not (0 < ro < rc):
        raise ValueError(f"LJ ro_A must be between 0 and rc_A (ro={ro}, rc={rc}).")

    return {"epsilon_eV": float(eps), "sigma_A": float(sig), "rc_A": float(rc), "ro_A": float(ro)}


def _calcMaxRc(atoms, margin=1e-3):
    import numpy as np

    a, b, c, alpha, betta, gamma = atoms.cell.cellpar()
    pbc = atoms.get_pbc()
    periodic_lengths = [L for L, is_p in zip((a, b, c), pbc) if is_p] #if not periodic, just use rc as it is

    if not periodic_lengths:
        return float('inf')

    L_min = min(periodic_lengths)
    return  0.4 * L_min    #L>2*rcut








