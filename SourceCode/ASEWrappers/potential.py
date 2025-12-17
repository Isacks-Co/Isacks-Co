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


import os
from asap3 import EMT
from asap3 import LennardJones as asap_LJ
class Potential:
    """Base class for interatomic potentials.

    This abstract interface defines the structure for all potentials
    that can produce ASE-compatible calculator objects.
    """

    def __init__(self):
        """Initialize the Potential base class."""
        self.pot_str = None

    def getASEPotentialCalculator(self):
        """Return an ASE calculator object.

        Returns:
            ase.Calculator | None: ASE-compatible calculator, or None if
            not implemented by the subclass.
        """
        return None

    def __str__(self):
        """Return the name of the potential.

        Returns:
            str: Human-readable potential name.
        """
        return self.pot_str


class LennardJonesPotential(Potential):
    """Lennard-Jones pair potential."""

    def __init__(self, atomic_numbers, sigmas, epsilons, rc=None):
        """Initialize Lennard-Jones parameters.

        Args:
            atomic_numbers (list[int]): List of atomic numbers.
            sigmas (list[float]): LJ sigma parameters (Å).
            epsilons (list[float]): LJ epsilon parameters (eV).
            rc (float, optional): Cutoff radius (Å). Defaults to `2.5 * max(sigmas)`.
        """
        super().__init__()
        self.pot_str = "Lennard Jones"

        self.atomic_nums = atomic_numbers
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.rc = rc if rc is not None else 2.5 * max(self.sigmas)

    def getASEPotentialCalculator(self):
        """Create an ASE Lennard-Jones calculator.

        Returns:
            ase.Calculator: ASAP LJ calculator instance.
        """
        return asap_LJ(
            self.atomic_nums,
            sigma=self.sigmas,
            epsilon=self.epsilons,
            rCut=self.rc
        )


class EMTPotential(Potential):
    """Effective Medium Theory (EMT) potential."""

    def __init__(self):
        """Initialize EMT potential."""
        super().__init__()
        self.pot_str = "EMT"

    def getASEPotentialCalculator(self):
        """Return the EMT ASE calculator.

        Returns:
            ase.calculators.emt.EMT: EMT calculator instance.
        """
        return EMT()


class MACEPotential(Potential):
    """MACE machine-learning interatomic potential."""

    def __init__(
        self,
        model_path=os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))
        ) + "/MACEModels/mace-mpa-0-medium.model"
    ):
        """Initialize the MACE potential.

        Args:
            model_path (str): Path to the .model file for MACE.
                Defaults to a packaged model in `MACEModels/`.
        """
        super().__init__()
        self.pot_str = "MACE"
        self.model_path = model_path

    def getASEPotentialCalculator(self):
        """Load and return the MACE ASE calculator.

        Includes internal environment variable settings and warning suppression.

        Returns:
            mace.calculators.MACECalculator: MACE calculator instance.
        """
        import os, warnings, torch

        os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

        warnings.filterwarnings(
            "ignore",
            message=(
                "Environment variable TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD detected, "
                "since the`weights_only` argument was not explicitly passed to "
                "`torch.load`, forcing weights_only=False."
            ),
        )
        warnings.filterwarnings(
            "ignore",
            message="cuequivariance or cuequivariance_torch is not available. "
                    "Cuequivariance acceleration will be disabled."
        )

        from mace.calculators import MACECalculator
        device = "cuda" if torch.cuda.is_available() else "cpu"

        return MACECalculator(
            model_paths=self.model_path,
            device=device,
            default_dtype="float64",
            head="default"
        )
