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
    def __init__(self, timestep):
        self.timestep = timestep
        self.atoms = None
        self.integrator_partial = None
        self.attachments = []

    @property
    def ensemble(self):
        return None

    def attach(self, attachment, interval):
        self.attachments.append((attachment, interval))

    def _addAttachments(self, integrator):
        for func, interval in self.attachments:
            integrator.attach(func, interval=interval)

    def run(self, atomic_structure, steps):
        with suppress_cpp_output():
            integrator_func = deepcopy(self.integrator_partial)(atomic_structure.getAtoms())

            self._addAttachments(integrator_func)

            integrator_func.run(steps)
        self.clearData()

    def clearData(self):
        self.atoms = None
        self.attachments = []

    def __str__(self):  # TODO Expand on this
        return self.ensemble


class VelocityVerletIntegrator(Integrator):
    def __init__(self, timestep):
        super().__init__(timestep=timestep)
        self.integrator_partial = partial(VelocityVerlet, timestep=self.timestep * fs)

    @property
    def ensemble(self):
        return "NVE"


class LangevinIntegrator(Integrator):
    def __init__(self, timestep, temperature_K, friction):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.friction = friction
        self.integrator_partial = partial(Langevin, timestep=self.timestep * fs, temperature_K=self.temperature_K,
                                          friction=self.friction / fs)

    @property
    def ensemble(self):
        return "NVT"


class IsotropicMTKNPTIntegrator(Integrator):
    def __init__(self, timestep, temperature_K, pressure, tdamp, pdamp):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.tdamp = tdamp
        self.pdamp = pdamp
        self.integrator_partial = partial(IsotropicMTKNPT, timestep=self.timestep * fs,
                                          temperature_K=self.temperature_K,
                                          pressure_au=self.pressure * GPa, tdamp=self.tdamp, pdamp=self.pdamp)

    @property
    def ensemble(self):
        return "NPT"


class BerendsenNPTIntegrator(Integrator):
    def __init__(self, timestep, temperature_K, pressure, compressibility):
        super().__init__(timestep=timestep)
        self.temperature_K = temperature_K
        self.pressure = pressure
        self.compressibility = compressibility
        self.integrator_partial = partial(NPTBerendsen, timestep=self.timestep * fs, temperature_K=self.temperature_K,
                                          pressure_au=self.pressure * GPa, compressibility_au=self.compressibility / GPa)

    @property
    def ensemble(self):
        return "NPT"
