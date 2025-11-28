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


def auToGPascal(pressure_au):
    """Convert pressure from eV/Å^3 to GPa."""
    return pressure_au * 1.602176634e2


def auToPascal(pressure_au):
    return auToGPascal(pressure_au) * 1e9


def pascalToAu(pressure_Pa):
    """Convert pressure from Pa to eV/Å^3."""
    return pressure_Pa * 6.241509074e-12


def GPascalToAu(pressure_GPa):
    return pressure_GPa * 6.241509074e-3


def specificHeatAuToSI(Cv_au):
    """Convert specific heat from eV/(amu·K) to J/(kg·K)."""
    return Cv_au * 9.64853321e7


def selfDiffusionCoeffAuToSI(D_au):
    return D_au * 1e-5


def atomicMassTokg(mass_au):
    return mass_au * 1.66053906660e-27


def kgToAtomicMass(mass_kg):
    return mass_kg / 1.66053906660e-27


def a1ToM1(length_a):
    return length_a * 1e-10


def mToA(length_m):
    return length_m * 1e10


def a2ToM2(area_a2):
    return area_a2 * 1e-20


def m2ToA2(area_m2):
    return area_m2 * 1e20


def a3ToM3(volume_a3):
    return volume_a3 * 1e-30


def m3ToA3(volume_m3):
    return volume_m3 * 1e30


def evToJ(energy_eV):
    return energy_eV * 1.602176634e-19


def jToEv(energy_J):
    return energy_J * 6.241509074e18


def hartreePerbohr3ToJPerM3(energy_hartreePerbohr3):
    return energy_hartreePerbohr3 * 2.19474631370215e13
