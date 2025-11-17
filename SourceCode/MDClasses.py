from ase.units import fs
from ase.io.trajectory import Trajectory
import numpy as np

from copy import copy

from ASEWrappers import DataTrajectory,Frame
from ASEWrappers import AtomicStructure
from Utils.equilibriumCondition import EquilibriumCondition
class MDBase:
    def __init__(self,settings):
        self.integrator = settings.integrator
        self.sample_data = None
        self.run_type = None

    def _storeFrame(self,atomic_strucuture: AtomicStructure,data_traj: DataTrajectory):
        if len(data_traj) == 0:
            time = 0
        else:
            time = data_traj[-1].time + self.integrator.timestep
      
        data = {label:  self._getAtomsData(atomic_strucuture,label,data_traj.initial_atoms) for label in  self.sample_data }

        frame = Frame(time,data)
 
        data_traj.append(frame)
  
    def _SaveASETrajectory(self,atomic_structure: AtomicStructure,interval = 1):
        traj = Trajectory(filename=f"{self.run_type}.traj", mode="w",atoms = atomic_structure.getAtoms())
        self.integrator.attach(traj.write,interval)

    @staticmethod
    def _getAtomsData(atomic_structure: AtomicStructure, name,initial_atomic_structure:AtomicStructure):
        """
        Help function for getting specific data from the ASE atoms object.
        """
     
        if atomic_structure.potential is None:
            raise RuntimeError("Atoms object has no potential")
        
        if name == "E_pot":
            E_pot = atomic_structure.potential_energy
            return E_pot

        if name == "E_kin":
            E_kin = atomic_structure.kinetic_energy
            return E_kin

        if name == "E_tot":
            E_tot = atomic_structure.total_energy
            return E_tot

        if name == "V":
            vol = atomic_structure.volume
            return vol
        
        if name == "T":
            T = atomic_structure.temperature
            return T

        if name == "MSD":
            MSD = atomic_structure.computeMSD(initial_atomic_structure)
            return MSD
        if name == "E_coh":
            E_coh = atomic_structure.cohesive_energy
            return E_coh
class EquilibriumRun(MDBase):
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Equil"
        self.sample_data = ["E_pot"]
    
    def run(self,atomic_structure: AtomicStructure ,num_steps,init_vel = False,store_traj = True):
        
        if init_vel:
            atomic_structure.setVelocitiesMB(self.integrator.temperature_K)
        
        if store_traj:
            self._SaveASETrajectory(atomic_structure)
        #self.integrator.attach(lambda: self._storeFrame(atomic_strucuture=atomic_structure,data_traj = data_traj),1)
        #self.integrator.attach(lambda: self._check_equilibrium(dataframe=data_traj),1)
   
        self.integrator.run(atomic_structure,num_steps)
        
       
        return atomic_structure
    
    def _check_equilibrium(self,dataframe):
        
        if len(dataframe) > 1000:
            if EquilibriumCondition.checkStable(dataframe[-1000:]):
                
                raise StopIteration("Equil reached")

class SampleRun(MDBase):
    def __init__(self, settings,sample_data = "all"):
        super().__init__(settings)
        self.run_type = "Sample"
        self.sample_data = ["T","E_tot","E_kin","E_pot","V","MSD"] if sample_data =="all" else sample_data
    def run(self,atomic_structure: AtomicStructure ,num_steps,store_traj = False):
        
        
        data_traj = DataTrajectory(atomic_structure)
        if store_traj:

            self._SaveASETrajectory(atomic_structure)

        self.integrator.attach(lambda: self._storeFrame(atomic_strucuture=atomic_structure,data_traj = data_traj),1)
        
        
        self.integrator.run(atomic_structure,num_steps) 
        # TODO Add equil check / Fail check
        return data_traj
    

class StrecthRun(MDBase): #TODO Finish this
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Stretch"
        
    def run(self,atomic_structure: AtomicStructure ): 
            
        strains = np.linspace(-0.005, 0.005, 5) # TODO Not hardcoded ? 
        cell0 = atomic_structure.cell
        stress0 = atomic_structure.stress
        hold_steps = 500 # TODO Not hardcoded ? 
        equil_atoms = copy(atomic_structure)
        calculator = atomic_structure.potential
        
        
        C = np.zeros((6, 6))
        for beta in range(6):
            
            # list for storing the average matrix of stresses for each strain.
            average_stress = []
            for e in strains:
                # Strain tensor in Voigt form
                stress_list = []
                eps = np.zeros((3, 3))
                if beta < 3:
                    eps[beta, beta] = e
                elif beta == 3:
                    eps[1, 2] = eps[2, 1] = e / 2.0
                elif beta == 4:
                    eps[0, 2] = eps[2, 0] = e / 2.0
                elif beta == 5:
                    eps[0, 1] = eps[1, 0] = e / 2.0

                # Apply the strain to the cell and perform the number of steps specified with hold_steps
                new_cell = np.dot(cell0, np.eye(3) + eps)
                atoms = copy(equil_atoms)
                #atoms.calc = calculator
                atoms.cell = new_cell
                self.integrator.attach(lambda : self.appendStress(atoms, stress_list, stress0),1)
                self.integrator.run(atoms,hold_steps)
                
                stacked = np.stack(stress_list)
                avg_matrix = stacked.mean(axis=0)
                average_stress.append(avg_matrix)
            sigmas = np.array(average_stress)
            for alpha in range(6):
                C[alpha, beta] = np.polyfit(strains, sigmas[:, alpha], 1)[0]
        np.save(f"cmatrix", C)
        return C
    

    def appendStress(self, atoms, stress_list, stress0):
        # Help function for appending stresses during _stretchCell runs
        stress = atoms.stress - stress0
    
        stress_list.append(stress)





