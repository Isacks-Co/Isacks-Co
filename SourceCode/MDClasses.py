from ase.units import fs
from ase.io.trajectory import Trajectory

from SourceCode.ASEWrappers import DataTrajectory,Frame
from SourceCode.ASEWrappers import AtomicStructure
from functools import partial

class MDBase:
    def __init__(self,settings):
        self.integrator = settings.integrator
        self.sample_data = None
        self.run_type = None

        
    def run(self,atomic_structure: AtomicStructure ,num_steps, store_traj, **kwargs):
        pass


    def _storeFrame(self,atomic_strucuture: AtomicStructure,data_traj: DataTrajectory):
        data = {label: self._getAtomsData(atomic_strucuture,label) for label in  self.sample_data }
        frame = Frame(data)
        data_traj.append(frame)

    def _SaveASETrajectory(self,atomic_structure: AtomicStructure,interval = 1):
        traj = Trajectory(filename=f"{atomic_structure.label}_{self.run_type}.traj", mode="w")
        self.integrator.attach(traj.write,interval)

    @staticmethod
    def _getAtomsData(atomic_structure: AtomicStructure, name):
        """
        Help function for getting specific data from the ASE atoms object.
        """
        if atomic_structure.potential is None:
            raise RuntimeError("Atoms object has no potential")
        
        if name == "E_pot":
            E_pot = atomic_structure.potential_energy()
            return E_pot

        if name == "E_kin":
            E_kin = atomic_structure.kinetic_energy()
            return E_kin

        if name == "E_tot":
            E_tot = atomic_structure.total_energy()
            return E_tot

        if name == "V":
            vol = atomic_structure.volume()
            return vol
        
        if name == "T":
            T = atomic_structure.temperature()
            return T

        if name == "F":
            F = atomic_structure.forces()
            return F
            
class EquilibriumRun(MDBase):
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Equil"
    
    def run(self,atomic_structure: AtomicStructure ,num_steps,store_traj = True):

        #TODO Should init_vel be here or is this "user responsibility" 

        if store_traj:
            self._SaveASETrajectory(atomic_structure)
        
        self.integrator.run(atomic_structure,num_steps) 
        # TODO Add equil check / Fail check
        return atomic_structure



class SampleRun(MDBase):
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Sample"
        
    def run(self,atomic_structure: AtomicStructure ,num_steps,store_traj = False):
        
        #TODO Should init_vel be here or is this "user responsibility" 
        data_traj = DataTrajectory(atomic_structure)
        if store_traj:
            self._SaveASETrajectory(atomic_structure)

        self.integrator.attach(self._storeFrame(atomic_strucuture=atomic_structure,data_traj = data_traj))
        
        
        self.integrator.run(atomic_structure,num_steps) 
        # TODO Add equil check / Fail check
        return data_traj
    

class StrecthRun(MDBase): #TODO Finish this
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Stretch"
        
    def run(self,atomic_structure: AtomicStructure ,num_steps,store_traj = True):
        pass
    

    

