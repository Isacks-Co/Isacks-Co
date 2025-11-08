from ase.units import fs

from DataContainer import DataTrajectory,Frame

class MDBase:
    def __init__(self,settings):

        self.timestep = settings.timestep * fs
        
        self.sample_data = None
        

        
    def run(self,atoms,num_steps, store_traj, **kwargs):
        pass


    def _storeFrame(self,atoms,data_traj: DataTrajectory):
        data = {label: self._getAtomsData(atoms,label) for label in  self.sample_data }
        frame = Frame(data)
        data_traj.append(frame)
    
    def _getAtomsData(atoms, name):
        """
        Help function for getting specific data from the ASE atoms object.
        """
        if atoms.calc is None:
            raise RuntimeError("Atoms object has no calculator")
        
        if name == "E_pot":
            E_pot = atoms.get_potential_energy()
            return E_pot

        if name == "E_kin":
            E_kin = atoms.get_kinetic_energy()
            return E_kin

        if name == "E_tot":
            E_tot = atoms.get_total_energy()
            return E_tot

        if name == "V":
            vol = atoms.get_volume()
            return vol
        
        if name == "T":
            T = atoms.get_temperature()
            return T

        if name == "F":
            F = atoms.get_forces()
            return F
            
class EquilibriumRun(MDBase):
    
    def run(self,atomic_structure ,num_steps, init_vel = True,store_traj = True):
        
        

        MaxwellBoltzmannDistribution(atoms, temperature_K=self.temperature_k,
                                     force_temp=True)  # Initialize velocity according to temperature_k
        Stationary(atoms) # Make sure center of mass has no linear momentum
        ZeroRotation(atoms) # Make sure center of mass has no angular momentum, might not be needed
        self.equilibriumRun(atoms=atoms) # TODO BREAKS TO EARLY

        log.info("MD run starts with: %i steps", self.steps)
        dyn = self.integrator(atoms=atoms)

        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms) ## currently have .. before

        dyn.attach(lambda: self.save_data(atoms,traj),
                   interval=self.interval)

        # Continue with the main MD run
        dyn.run(self.steps)  # RUN




