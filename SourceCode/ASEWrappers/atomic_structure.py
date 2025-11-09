
import numpy as np
from ase import Atoms
from ase.io import read
import hashlib
from SourceCode.ASEWrappers.potential import Potential
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary, ZeroRotation


class AtomicStructure:
    """
    Wrapper around ASE Atoms
    """

    def __init__(self, atoms: Atoms, special_label = None):
        self._atoms = atoms.deepcopy() 

        self._label = self._generateHashLabel(special_label)

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
    def potential_energy(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get energy: no potential assigned")
        return self._atoms.get_potential_energy()

    @property
    def kinetic_energy(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get energy: no potential assigned")
        return self._atoms.get_kinetic_energy()
    
    @property
    def total_energy(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get energy: no potential assigned")
        return self._atoms.get_total_energy()
    
    @property
    def temperature(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get energy: no potential assigned")
        return self._atoms.get_temperature()
    
    @property
    def volume(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get energy: no potential assigned")
        return self._atoms.get_volume()
    
    @property
    def cell(self):
        return self._atoms.cell.array.copy()

    @cell.setter
    def cell(self, new_cell):
        self._atoms.set_cell(new_cell, scale_atoms=True)

    @property
    def symbols(self):
        return self._atoms.get_chemical_symbols()
    



    def setVelocitiesMB(self,temperature_K):
        MaxwellBoltzmannDistribution(self._atoms, temperature_K=self.temperature_K,
                                     force_temp=True)  # Initialize velocity according to temperature_k
        Stationary(self._atoms) # Make sure center of mass has no linear momentum
        ZeroRotation(self._atoms) # Make sure center of mass has no angular momentum, might not be needed

    def getAtoms(self):
        return self._atoms
    
    


    @property
    def label(self):
        return self._label 
    
    def _generateHashLabel(self,special_label):
        

        formula = self._atoms.get_chemical_formula()
        data = (self._atoms.numbers.tobytes() + self._atoms.positions.tobytes())
        short_hash = hashlib.md5(data).hexdigest()[:6]  # first 6 chars
        if special_label == None:
            label = f"{formula}_{short_hash}"
        else:
            label = f"{formula}_{special_label}_{short_hash}"
        return label
