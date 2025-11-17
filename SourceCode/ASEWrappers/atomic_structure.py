import numpy as np
from ase import Atoms
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary, ZeroRotation

import hashlib

from .potential import Potential


class AtomicStructure:
    """
    Wrapper around ASE Atoms.
    Can be created from eithre a complete atoms object, including potential or directly 
    from a file in combination with a potential object.
    """

    def __init__(self, atoms: Atoms, special_label = None,label = None):
        self._atoms = atoms.copy()
        if atoms.calc != None:
            self._atoms.calc = atoms.calc
        else:
            raise ValueError("Atoms object needs to have an attached calculator")
        self._label = self._generateHashLabel(special_label) if label == None else label

    @classmethod
    def fromFile(cls,path, pbc = True,supercells = [1,1,1], potential:Potential = None):
        atoms = read(path) * supercells
        atoms.calc = potential.getASEPotentialCalculator()
        atoms.pbc = pbc
        return cls(atoms)

    #Standard dunder methods

    def __len__(self):
        return len(self._atoms)
    
    def __str__(self):
        return f"{self._atoms}"

    def __copy__(self):
        return AtomicStructure(self._atoms,label=self.label)
    
    # Properties. 

    @property
    def label(self):
        return self._label 
    
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
    def stress(self):
        if self._atoms.calc is None:
            raise RuntimeError("Cannot get forces: no potential assigned")
        return self._atoms.get_stress(voigt=True)

    @property
    def masses(self):
        return self._atoms.get_masses()
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

    @property
    def cohesive_energy(self):
        
        """
        Calculate the cohesive energy per atom. 
        This is done by computing the difference in energy between the separate atoms and the bulk structure 
        
        Unit: ev/atom
        """
        
        number = self.__len__()
        e_atoms = 0

        for symbol in self._atoms.get_chemical_symbols():
            atom = Atoms(symbol, positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=False)

            atom.calc = self._atoms.calc
            e_atoms += atom.get_potential_energy()
        
        e_bulk =self.potential_energy
        
        e_coh = (e_atoms-e_bulk) / number

        return e_coh

    
    @property
    def internal_pressure(self):
        """
        Compute internal pressure using atomic units internally.
        Instantaneous: P = (1/3V) [ 2 N E_kin + sum_i r_i · f_i ]
        Returns a mean over the instantaneous frames in the trajectory
        Unit: ev/Å^3
        """

        N = len(self)

        
        e_kin_eV = self.kinetic_energy
        V_A3 = self.volume
        forces_eVA = self.forces
        positions_A = self.positions
        sum_rf = np.sum(forces_eVA * positions_A)
        P_eVA3 = (1.0 / (3.0 * V_A3)) * (2.0 * N * e_kin_eV + sum_rf)
        
        return P_eVA3


    # Help functions 

    def setVelocitiesMB(self,temperature_K):
        """
        Sets the velocities of the structure at a given temperature accordint to a 
        Maxwell-Boltzmann distribution. Also makes sure that the result
        has no drift.

        In:
            (float) temperature_K: wanted temperature of the structure in Kelvin
        """

        MaxwellBoltzmannDistribution(self._atoms, temperature_K=temperature_K,
                                     force_temp=True)  # Initialize velocity according to temperature_k
        Stationary(self._atoms) # Make sure center of mass has no linear momentum
        ZeroRotation(self._atoms) # Make sure center of mass has no angular momentum, might not be needed

    def computeMSD(self,orig_struct):
        """
        Computes MSD in relation to some reference structure
        In:
            (AtomicStructure) : reference structure

        Out: 
            (int) : Mean square displacement over all atoms
        """
        if not isinstance(orig_struct,AtomicStructure):
            raise TypeError("orig_struct needs to be type AtomicStructure")
        r_0 = orig_struct.positions
        r_n = self.positions
        return np.mean((r_0 - r_n) ** 2)
    
    def _generateHashLabel(self,special_label):
        """
        Generates a semi-random hashlabel for the atomic structure.
        In:
            (str) special label : custom label added in the middle of the hashlabel
        Out:
            (str) hashlabel on the form (chemical formula)_(special label)_(semi-random hash)
        """
        formula = self._atoms.get_chemical_formula()
        data = (self._atoms.numbers.tobytes() + self._atoms.positions.tobytes() + np.random.bytes(4))
        short_hash = hashlib.md5(data).hexdigest()[:6]  
        if special_label == None:
            label = f"{formula}_{short_hash}"
        else:
            label = f"{formula}_{special_label}_{short_hash}"
        return label
    
    def getAtoms(self): # help function that should be used carefully
        return self._atoms