import numpy as np
from ase import Atoms
from ase.units import fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase.md import MDLogger
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
import functools



atoms = FaceCenteredCubic(size=(1,1,1), symbol="Cu", pbc=True)


class MDbase:
    """
        basic MD class
        from preprocessing class we should get an atom-object with initialized settings.
        When initialized it represents a MD-simulation with prefilled settings that can be used for multiple runs with 
        different atomic structures.
    """

    def __init__(self, ts_fs, ts, interval,
                  integrator_str, r_path,
                  T_k, friction, potential_str,att_list,
                  pressure,comp = 10**(-11)):
        """
        In:
            ts_fs : timestesp (femto) FLOAT
            ts : timesteps INT
            inteval: dead variable will be removing soon.
            r_path : output file path STR
            T_k : Temperature in Kelvin for Langevin simualtion
            friction: Friction for Langevin simulation
            potential_str: potential to be used in simulation.
            att_list: list of attachments
            
        """
        self.dt = float(ts_fs * fs)
        self.steps = int(ts)
        self.interval = int(interval) if interval else 10
        self.raw_path = r_path
        self.T_k = float(T_k)
        self.friction = float(friction)
        self.pressure = float(pressure)
        self.compressibility = comp
        self.potential = self.get_potential(potential_str)
        self.integrator = self.get_integrator(integrator_str)
        self.attachments = self.get_attachment(att_list)

        
    
    def get_potential(self, pot_str):
        """
        In:
            String: pot_str
        Out:
            Potential_function: potential
        """
        if pot_str in ["emt", "EMT"]:
            from asap3 import EMT as asap_EMT
            return asap_EMT
        
        elif pot_str in ["LJ", "lj", "lennardjones", "LennardJones", "Lennard_Jones"]:
            from ase.calculators.lj import LennardJones
            return LennardJones



    def get_integrator(self, str):
        if str in ["verlet", "Verlet"]:
            from ase.md.verlet import VelocityVerlet #för NVE
            return functools.partial(VelocityVerlet, timestep = self.dt )
        
        elif str in ["langevin", "Langevin"]:
            from ase.md.langevin import Langevin #för NVT
            return functools.partial(Langevin, timestep=self.dt,  temperature_K=self.T_k,  friction = self.friction)
        

        elif str in ["Berendsen", "NPT"]:
            from ase.md.nptberendsen import NPTBerendsen
            return functools.partial(NPTBerendsen, timestep=self.dt,  temperature_K=self.T_k,  
                                     pressure_au = self.pressure,compressibility_au = self.compressibility)

        
        else:
            raise ValueError("Invalid integrator")

    def get_attachment(self, attachments):
        pos_attachments = {'energy':self.print_energy, 
                           "momenta":self.print_momentum,
                           "center_of_mass":self.print_center_of_mass}
    
        return [pos_attachments[a] for a in attachments]



    def run_MD(self, atoms):
        """
        In: 
            Atoms: ase Atoms object representing the crystal structure

        Runs a MD simulation with the setting specified in __init__
        Depending on attachments will possibly print some data.
        Will always save a trajectory and log file.        
        """

        atoms.calc = self.potential()

        MaxwellBoltzmannDistribution(atoms, temperature_K=self.T_k,force_temp=True) # Initialize velocity according to T_k
        
        dyn = self.integrator(atoms=atoms) 
        
        material_name = str(atoms.symbols)
        print("MATERIALNAMN: ", material_name)

        #save traj
        traj = Trajectory(filename= f"{self.raw_path}.traj", mode="w", atoms=atoms)

        for a in self.attachments: 
            dyn.attach(functools.partial(a, atoms = atoms),interval = self.interval) # Attach the different functions for printing
        
        dyn.attach(traj.write, interval=self.interval)

        logger = MDLogger(dyn, atoms=atoms, logfile=f"{self.raw_path}.log", 
                          header=True, peratom=True, mode='a') #Create a logger for writing data
        dyn.attach(logger, interval=self.interval) # Attach logger

        dyn.run(self.steps) #RUN
    
    def print_energy( self,atoms):
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        etot = epot + ekin
        T = float(atoms.get_temperature())

        print(f"E_pot/atom={epot:.5f}  E_kin/atom={ekin:.5f}  E_tot/atom={etot:.5f}  T={T:.1f} K")

    def print_momentum( self,atoms):
        momenta = atoms.get_momenta()
        T = float(atoms.get_temperature())
        print(f"momenta: {momenta}  T={T:.1f} K")

    
    def print_center_of_mass(self,atoms):
        momenta = atoms.get_center_of_mass()
        print(f"Center of mass: {momenta}")


    def print_lattice_constants(self,atoms):
        print("Lattice: ",atoms.cell.cellpar())
        
       
def vizualise_traj(traj_file):
    traj = Trajectory(traj_file)
    view(traj)   



if __name__ == "__main__":
    md = MDbase(
        ts_fs=0.5,
        ts=1000,
        interval=10,       # print var 10:e steg
        integrator_str="Berendsen", #används inte just nu
        r_path="data",
        T_k=300.0,
        friction=0.01, #ingen aning vad som är rimligt
        potential_str="EMT",
        att_list= [],
        pressure= 10e+6
    )
    md.run_MD(atoms)
    vizualise_traj("data.traj")
