# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


LJ_DB = {
    "ar": {"epsilon_eV": 0.0103, "sigma": 3.405},
    "kr": {"epsilon_eV": 0.01475, "sigma": 3.650},
    "xe": {"epsilon_eV": 0.01903, "sigma": 4.100},
    "ne": {"epsilon_eV": 0.00307, "sigma": 2.789},
    "cu": {"epsilon_eV": 0.40, "sigma": 2.28},
    "fe": {"epsilon_eV": 0.2007, "sigma": 2.4193},
    "ni": {"epsilon_eV": 0.1729, "sigma": 1.5808},
    "al": {"epsilon_eV": 0.10, "sigma": 2.61},
    "au": {"epsilon_eV": 0.23, "sigma": 2.57},
    "pt": {"epsilon_eV": 0.28, "sigma": 2.47},
    # These are for tests
    "cr": {"epsilon_eV": 0.67322, "sigma": 2.2813},
    "po": {"epsilon_eV": 0.2, "sigma": 3.6},
}


def LJParams(material: str, *, epsilon_eV=None, sigma_A=None, rc_A=None,
             ro_A=None):
    """
    Help function to set sigma, epsilon and rCutoff for the given materials in the LJ_DB dictionary.
    Args:
        material (str): The material to set sigma, epsilon and rCutoff
        epsilon_eV (float, optional): The epsilon eVf rom the potential well
        sigma_A (float, optional): The sigma A from the potential well
        rc_A (float, optional): The cutoff radius
        ro_A (float, optional): The smooth of radius
    Returns:
        dict: The sigma, epsilon and rCutoff values

    """
    if material != "":
        base = LJ_DB.get(material.lower())
        if not base:
            raise ValueError(f"No LJ-parametrar registered for '{material}'."
                             f"Change to any of these: "
                             f"{', '.join(LJ_DB.keys())}"
                             )
    else:
        base = {}

    eps = epsilon_eV if epsilon_eV is not None else base.get("epsilon_eV")
    sig = sigma_A if sigma_A is not None else base.get("sigma")
    if eps is None or sig is None:
        raise ValueError("Lennard-Jones needs epsilon_eV and sigma.")

    rc = rc_A if rc_A is not None else 2.5 * float(sig)
    ro = ro_A if ro_A is not None else 0.9 * float(rc)

    if not (0 < ro < rc):
        raise ValueError(
            f"LJ ro_A must be between 0 and rc_A (ro={ro}, rc={rc}).")  # Doesnt really make sense as error message since user have no control of these values
    return {"epsilon_eV": float(eps), "sigma_A": float(sig), "rc_A": float(rc), "ro_A": float(ro)}


def calcMaxRc(atoms, margin=1e-3):
    """
    Function to calculate the maximal rcutoff allowed.
    Args:
        atoms : the atomobject to calculate the radius from

    Returns:
        float: The maximum rcutoff allowed
    """
    a, b, c, alpha, betta, gamma = atoms.cell.cellpar()
    pbc = atoms.get_pbc()
    periodic_lengths = [L for L, is_p in zip((a, b, c), pbc) if is_p]  # if not periodic, just use rc as it is

    if not periodic_lengths:
        return float('inf')

    L_min = min(periodic_lengths)
    return 0.4 * L_min  # L>2*rcut
