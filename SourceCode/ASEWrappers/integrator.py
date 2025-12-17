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


from asap3.md.langevin import Langevin
from asap3.md.nose_hoover_chain import IsotropicMTKNPT
from asap3.md.nptberendsen import NPTBerendsen
from asap3.md.verlet import VelocityVerlet
from ase.units import fs, GPa
from copy import deepcopy
from functools import partial

from .contextManager import suppress_cpp_output

"""
Wrapper of the ASE integrators. 
The purpose is to make the code more readable 
by wrapping up some of the complexities ASEs 
interface results in. 
"""

class Integrator:
    """Base class for all molecular dynamics integrators.

    This class defines the common interface for ASE-based integrators,
    including timestep handling, attachment of callback functions, and
    running integration loops.
    """

    def __init__(self, timestep):
        """Initialize the integrator.

        Args:
            timestep (float): Integration timestep (in femtoseconds).
        """
        self.timestep = timestep
        self.atoms = None
        self.integrator_partial = None
        self.attachments = []

    @property
    def ensemble(self):
        """str: Thermodynamic ensemble associated with the integrator."""
        return None

    def attach(self, attachment, interval):
        """Attach a callback function to the integrator.

        Args:
            attachment (callable): Function executed every `interval` steps.
            interval (int): Number of steps between calls.
        """
        self.attachments.append((attachment, interval))

    def _addAttachments(self, integrator):
        """Register all stored attachments with the ASE integrator.

        Args:
            integrator: ASE integrator instance that supports attach().
        """
        for func, interval in self.attachments:
            integrator.attach(func, interval=interval)

    def run(self, atomic_structure, steps):
        """Run the integrator on the provided atomic structure.

        Args:
            atomic_structure (AtomicStructure): Structure to integrate.
            steps (int): Number of MD steps to perform.
        """
        with suppress_cpp_output():
            integrator_func = deepcopy(self.integrator_partial)(
                atomic_structure.getAtoms()
            )

            self._addAttachments(integrator_func)
            integrator_func.run(steps)

        self.clearData()

    def clearData(self):
        """Reset stored atoms and clear attachment list."""
        self.atoms = None
        self.attachments = []

    def __str__(self):
        """Return the ensemble name.

        Returns:
            str: Ensemble string (e.g., "NVE", "NVT", "NPT").
        """
        return self.ensemble


class VelocityVerletIntegrator(Integrator):
    """Velocity-Verlet molecular dynamics integrator (NVE ensemble)."""

    def __init__(self, timestep):
        """Initialize Velocity Verlet integrator.

        Args:
            timestep (float): Integration timestep (fs).
        """
        super().__init__(timestep=timestep)
        self.integrator_partial = partial(VelocityVerlet, timestep=self.timestep * fs)

    @property
    def ensemble(self):
        """str: Thermodynamic ensemble (NVE)."""
        return "NVE"


class LangevinIntegrator(Integrator):
    """Langevin molecular dynamics integrator (NVT ensemble)."""

    def __init__(self, timestep, temperature_K, friction):
        """Initialize Langevin integrator.

        Args:
            timestep (float): Integration timestep (fs).
            temperature_K (float): Target temperature (K).
            friction (float): Friction coefficient (1/fs).
        """
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.friction = friction
        self.integrator_partial = partial(
            Langevin,
            timestep=self.timestep * fs,
            temperature_K=self.temperature_K,
            friction=self.friction / fs
        )

    @property
    def ensemble(self):
        """str: Thermodynamic ensemble (NVT)."""
        return "NVT"


class IsotropicMTKNPTIntegrator(Integrator):
    """MTK NPT integrator with isotropic cell fluctuations (NPT ensemble)."""

    def __init__(self, timestep, temperature_K, pressure, tdamp, pdamp):
        """Initialize MTK NPT integrator.

        Args:
            timestep (float): Integration timestep (fs).
            temperature_K (float): Target temperature (K).
            pressure (float): Target pressure (GPa).
            tdamp (float): Temperature damping time (fs).
            pdamp (float): Pressure damping time (fs).
        """
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.tdamp = tdamp
        self.pdamp = pdamp
        self.integrator_partial = partial(
            IsotropicMTKNPT,
            timestep=self.timestep * fs,
            temperature_K=self.temperature_K,
            pressure_au=self.pressure * GPa,
            tdamp=self.tdamp,
            pdamp=self.pdamp
        )

    @property
    def ensemble(self):
        """str: Thermodynamic ensemble (NPT)."""
        return "NPT"


class BerendsenNPTIntegrator(Integrator):
    """Berendsen NPT integrator (NPT ensemble)."""

    def __init__(self, timestep, temperature_K, pressure, compressibility):
        """Initialize Berendsen NPT integrator.

        Args:
            timestep (float): Integration timestep (fs).
            temperature_K (float): Target temperature (K).
            pressure (float): Target pressure (GPa).
            compressibility (float): Compressibility (1/GPa).
        """
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.compressibility = compressibility
        self.integrator_partial = partial(
            NPTBerendsen,
            timestep=self.timestep * fs,
            temperature_K=self.temperature_K,
            pressure_au=self.pressure * GPa,
            compressibility_au=self.compressibility / GPa
        )

    @property
    def ensemble(self):
        """str: Thermodynamic ensemble (NPT)."""
        return "NPT"
