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

from ase.build import bulk
from ase.visualize import view

from ASEWrappers.potential import MACEPotential

logger = logging.getLogger(__name__)


def getChemicalPotential(element: str):
    """
    Function to calculate the chemical potential of the given element, by using the approximated value from E_pot/N
    Args:
        element (str): The element from the periodic table
    """
    material_bulk = []
    try:
        material_bulk = bulk(element)
    except Exception:
        logger.error(f"No such element is accepted by ase.build.bulk")

    material_bulk = material_bulk.repeat(3)
    view(material_bulk)
    material_bulk.calc = MACEPotential().getASEPotentialCalculator()

    print(f"This is Chemical Potential for {element} : {material_bulk.get_potential_energy() / len(material_bulk)}")

    return material_bulk.get_potential_energy() / len(material_bulk)


if __name__ == "__main__":
    getChemicalPotential("Ta")
