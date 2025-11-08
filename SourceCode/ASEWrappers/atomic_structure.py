
import numpy as np
from ase import Atoms
from ase.io import read

from potential import Potential

class AtomicStructure:
    """
    Wrapper around ASE Atoms
    """

    def __init__(self, atoms: Atoms):
        self._atoms = atoms.copy() 

    @classmethod
    def fromFile(cls,path, pbc = True,supercells = [1,1,1], potential:Potential = None):

        atoms = read(path) * supercells
        atoms.calc = potential.getASEPotentialCalculator()
        atoms.pbc = pbc
        return cls(atoms)

    @property
    def potential(self):
        return self._atoms.calc
    
    @potential.setter
    def potential(self,new_pot):
        self._atoms.calc = new_pot
    
    @property
    def positions(self):
        return self._atoms.get_positions()

    @positions.setter
    def positions(self, new_positions):
        self._atoms.set_positions(new_positions)

    @property
    def velocities(self):
        vel = self._atoms.get_velocities()
        return vel.copy() if vel is not None else None

    @velocities.setter
    def velocities(self, new_vel):
        self._atoms.set_velocities(new_vel)

    @property
    def forces(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get forces: no potential assigned")
        return self._atoms.get_forces()

    @property
    def energy(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get energy: no potential assigned")
        return self._atoms.get_potential_energy()

    @property
    def cell(self):
        return self._atoms.cell.array.copy()

    @cell.setter
    def cell(self, new_cell):
        self._atoms.set_cell(new_cell, scale_atoms=True)

    @property
    def symbols(self):
        return self._atoms.get_chemical_symbols()

    def get_atoms(self):
        return self._atoms