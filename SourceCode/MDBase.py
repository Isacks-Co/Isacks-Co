import functools
import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import fs
from ase.visualize import view
from SourceCode.logger import logger_setup

log = logger_setup()
class MDBase:
    """
        basic MD class
        from preprocessing class we should get an atom-object with initialized settings.
        When initialized it represents a MD-simulation with prefilled settings that can be used for multiple runs with 
        different atomic structures.
    """

    def __init__(self, timestep_fs: float = 2, number_of_steps: int = 200, interval: int = 10,
                 integrator_str: str = "Verlet", output_file: str = "data",
                 temperature_k: float = 293, friction: float = 0.01, potential_str: str = "EMT",
                 att_list: list = ["energy"],
                 pressure: float = 10e+6, compressibility: float = 10e-11,  equil_steps: int = 2000):
        """
        In:
            timestep_fs : timestesp (femto) FLOAT
            number_of_steps : timesteps INT
            inteval: dead variable will be removing soon.
            output_file : output file path STR
            temperature_k : Temperature in Kelvin for Langevin simualtion
            friction: Friction for Langevin simulation
            potential_str: potential to be used in simulation.
            att_list: list of attachments
            
        """
        self.timestep = float(timestep_fs * fs)
        self.steps = int(number_of_steps)
        self.interval = int(interval) if interval else 10
        self.output_file = output_file
        self.temperature_k = float(temperature_k)
        self.friction = float(friction) / fs
        self.pressure = float(pressure)
        self.compressibility = compressibility
        self.potential = self.getPotential(potential_str)
        self.integrator = self.getIntegrator(integrator_str)
        self.attachments = self.getAttachment(att_list)
        self.equilibrium_steps = equil_steps
        self.temp_history = []
        self.hits = 0
        self.ensemble = integrator_str

        log.debug(
        "MDBase init: dt(fs)=%s steps=%s interval=%s T=%sK friction=%s pot=%s integrator=%s out=%s",
        timestep_fs, number_of_steps, self.interval, self.temperature_k, self.friction,
        potential_str, integrator_str, self.output_file
        )


    @classmethod
    def initNVE(cls, temperature: float,  pot_str:str, timestep:float,
                steps:int, interval:int,output_file: str, equilibrium_steps:int):

        return cls(temperature_k = temperature, integrator_str = "NVE", potential_str = pot_str ,
                timestep_fs = timestep, number_of_steps = steps, interval = interval, equil_steps = equilibrium_steps, output_file = output_file)

    @classmethod
    def initNVT(cls, temperature: float, friction: float,  pot_str:str, timestep:float,
                steps:int, interval:int, equilibrium_steps:int, output_file:str):
        return cls(temperature_k = temperature, friction = friction,  integrator_str = "NVT", potential_str = pot_str,
                    timestep_fs = timestep, number_of_steps = steps, interval = interval, equil_steps = equilibrium_steps, output_file=output_file )

    @classmethod
    def initNPT(cls, temperature: float, timestep: float,
                steps: int, interval: int, pressure_Pa: float, compressibility: float, pot_str: str, equilibrium_steps: int, output_file: str):
        return cls(temperature_k=temperature, pressure=pressure_Pa, compressibility=compressibility,
                   integrator_str="NPT", potential_str=pot_str, timestep_fs=timestep,
                   number_of_steps=steps, interval=interval,equil_steps=equilibrium_steps, output_file=output_file)


    def pascalToAu(self, pressure_Pa):
        pressure_au = pressure_Pa * 6.2415e-12
        return pressure_au



    def getPotential(self, potential: str):
        """
        In:
            String: potential
        Out:
            Potential_function: potential
        """
        potential_lower = potential.lower()
        if potential_lower in ["emt"]:
            from asap3 import EMT as asap_EMT
            log.info("Potential: EMT")
            return asap_EMT

        elif potential_lower in ["lj", "lennardjones", "lennard_jones"]:
            from ase.calculators.lj import LennardJones
            log.info("Potential: Lennard Jones")
            return LennardJones
        else:
            log.error("Invalid potential function: %s", potential)
            raise ValueError(f"Invalid potential function: {potential}")

    def getIntegrator(self, integrator: str):
        integrator_lower = integrator.lower()
        if integrator_lower in ["verlet", "nve"]:
            from asap3.md.verlet import VelocityVerlet  # för NVE
            log.info("Integrator: Verlet")
            return functools.partial(VelocityVerlet, timestep=self.timestep)

        elif integrator_lower in ["langevin", "nvt"]:
            from asap3.md.langevin import Langevin  # för NVT
            log.info("Integrator: Langevin")
            return functools.partial(Langevin, timestep=self.timestep, temperature_K=self.temperature_k, friction=self.friction)

        elif integrator_lower in ["berendsen", "npt"]:
            from asap3.md.nptberendsen import NPTBerendsen
            log.info("Integrator: Berendsen")
            return functools.partial(NPTBerendsen, timestep=self.timestep, temperature_K=self.temperature_k,
                                     pressure_au=self.pressure, compressibility_au=self.compressibility)

        else:
            log.error("Invalid Integrator function: %s", integrator) ##
            raise ValueError(f"Invalid integrator: {integrator}")


    def getAttachment(self, attachments):
        pos_attachments = {'energy': self.printEnergy,
                           "momenta": self.printMomentum,
                           "center_of_mass": self.printCenterOfMass,
                           "lattice":self.printLatticeConstants }
        
        for a in attachments:
            if a not in pos_attachments.keys():
                raise ValueError(f"Invalid attachment: {a}")


        return [pos_attachments[a] for a in attachments]

    def equilibriumRun(self, atoms):
                     
        #NVT until equilibrium is reached
        from asap3.md.langevin import Langevin
        dyn_eq = Langevin(atoms,
                          timestep=self.timestep,
                          temperature_K=self.temperature_k,
                          friction=self.friction)
        

        #traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms)
        #dyn_eq.attach(traj.write, interval=self.interval)
        log.info(f"Starting equilibrium run with NVT Ensemble to reach desired temperature of {self.temperature_k} K")

        dyn_eq.run(int(self.equilibrium_steps))
        current_T = atoms.get_temperature()
        log.info(f"Systems temperature is {round(current_T,2)} K after {self.equilibrium_steps} steps")
        


    def runMD(self, atoms):
        """
        In: 
            Atoms: ase Atoms object representing the crystal structure

        Runs a MD simulation with the setting specified in __init__
        Depending on attachments will possibly print some data.
        Will always save a trajectory and log file.        
        """

        atoms.calc = self.potential()

       
        MaxwellBoltzmannDistribution(atoms, temperature_K=self.temperature_k,
                                     force_temp=True)  # Initialize velocity according to temperature_k

        self.equilibriumRun(atoms=atoms)
        log.info("MD run starts with: %i steps", self.steps)

        dyn = self.integrator(atoms=atoms)

        #material_name = str(atoms.symbols)
        #print("MATERIALNAMN: ", material_name)

        # save traj
        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms)

        for a in self.attachments:
            dyn.attach(functools.partial(a, atoms=atoms),
                       interval=self.interval)  # Attach the different functions for printing

        dyn.attach(traj.write, interval=self.interval)

        logger = MDLogger(dyn, atoms=atoms, logfile=f"{self.output_file}.log",
                          header=True, peratom=True, mode='a')  # Create a logger for writing data
        dyn.attach(logger, interval=self.interval)  # Attach logger
        dyn.attach(lambda: self.failSafe(atoms), interval=self.interval)
        dyn.run(self.steps)  # RUN

    def printEnergy(self, atoms):
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        etot = (epot + ekin) 
        T = float(atoms.get_temperature())

        print(f"E_pot/atom={epot:.5f}  E_kin/atom={ekin:.5f}  E_tot/atom={etot:.5f}  T={T:.1f} K")

    def printMomentum(self, atoms):
        momenta = atoms.get_momenta()
        T = float(atoms.get_temperature())
        print(f"momenta: {momenta}  T={T:.1f} K")

    def printCenterOfMass(self, atoms):
        momenta = atoms.get_center_of_mass()
        print(f"Center of mass: {momenta}")

    def printLatticeConstants(self, atoms):
        print("Lattice: ", atoms.cell.cellpar())

    def visualizeTraj(self):
        traj = Trajectory("data.traj")
        view(traj)

    def failSafe(self, atoms):
        """
        Checks if temperature diverges continously in one direction,
        returns an error if that's the case, uses a 'window' number of last runs to calculate mean average
        hits determines how many flucuations before exiting
        """
        window = 20
        if (self.ensemble != "NVE"):
            self.temp_history.append(atoms.get_temperature())
        else:
            self.temp_history.append((atoms.get_potential_energy()+atoms.get_kinetic_energy())/len(atoms))
        if len(self.temp_history) > 1:
            mean = np.mean(self.temp_history[-window:])
            std  = np.std(self.temp_history[-window:])
            # TODO Changed for problematic effects, needs to be looked at
            # if abs(self.temp_history[0] - mean) > 2 * std:
            if std > self.temp_history[0]:
                if (self.hits == 2):
                    if (self.ensemble != "NVE"):
                        raise RuntimeWarning("Run canceled because simulation is not stable. Temperature change is greater than 2 standard deviations.")
                    else:
                        raise RuntimeWarning("Run canceled because simulation is not stable. Total energy change is greater than 2 standard deviations.")
                else:
                    self.hits += 1
                                