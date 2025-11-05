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

import logging
from LJRegistry import LJParams, calcMaxRc

log = logging.getLogger(__name__)


class Potential:
    
    def getPotential(self, potential: str):
        """
        In:
            String: potential
        Out:
            Potential_function: potential
        """
        potential_lower = potential.lower()
        if potential_lower in ["emt"]:
            log.info("Potential: EMT")
            return self.setupEMT

        elif potential_lower in ["lj", "lennardjones", "lennard_jones"]:
            log.info("Potential: Lennard Jones")
            return self.setupLJCalculator

        elif potential_lower in ["mace", "MACE"]:
            return self.setUpMACE

        else:
            log.error("Invalid potential function: %s", potential)
            raise ValueError(f"Invalid potential function: {potential}")

    def setUpMACE(self, atoms):
        import os, warnings
        os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
        warnings.filterwarnings(
            "ignore",
            message="Environment variable TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD detected",
        )

        from mace.calculators import MACECalculator
        import torch
        model_path = "MACEModels/mace-mpa-0-medium.model" #just to try
        device = "cuda" if torch.cuda.is_available() else "cpu"
        return MACECalculator(model_paths= model_path, device=device, default_dtype="float64",head="default")

    def setupEMT(self, atoms):
        from asap3 import EMT as asap_EMT
        # atoms not used, but lets us skip the lambda solution
        return asap_EMT()
        # return lambda atoms: asap_EMT()

    def setupLJCalculator(self, atoms):
        symbols = atoms.get_chemical_symbols()
        uniq = sorted(set(symbols))
        if len(uniq) != 1:
            raise ValueError(
                f"ASE LennardJones supports a single atom type only; found {uniq}. "
            )

        material_key = uniq[0].lower()
        params = LJParams(material=material_key)
        atomic_number = [(atoms.get_atomic_numbers()[0])]
        eps = params["epsilon_eV"]
        sig = params["sigma_A"]
        ro = params["ro_A"]
        rc = params["rc_A"]
        rc_max = calcMaxRc(atoms)

        if rc > rc_max:
            log.warning("The rCut is larger than the cell size, will use cell size to derive new value for rCut")
            rc = rc_max
            ro = 0.9 * float(rc)

        log.debug("Using  rc =  ", rc)
        try:
            from asap3 import LennardJones as asap_LJ
            calc_asap = asap_LJ(
                atomic_number,
                epsilon=[eps],
                sigma=[sig],
                rCut=rc,
                modified=True
            )

            # atoms.calc = calc_asap
            # _ = atoms.get_potential_energy()
            log.info("Using asap3 LJ | element=%s (Z=%s) | ε=%.4g eV | σ=%.4g Å | rc=%.4g Å ",
                     material_key, atomic_number[0], eps, sig, rc)
            return calc_asap


        except Exception as e:
            from ase.calculators.lj import LennardJones as ase_LJ
            calc_ase = ase_LJ(
                epsilon=eps,
                sigma=sig,
                rc=rc,
                ro=ro
            )
            log.warning(
                f"Falling back to ASE LJ | Reason: {e}"
            )
            return calc_ase
