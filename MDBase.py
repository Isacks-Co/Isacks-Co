import functools
import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import fs
from ase.visualize import view


class MDBase:
    """
        basic MD class
        from preprocessing class we should get an atom-object with initialized settings.
        When initialized it represents a MD-simulation with prefilled settings that can be used for multiple runs with 
        different atomic structures.
    """

    def __init__(self, timestep_fs: int = 10, number_of_steps: int = 200, interval: int = 10,
                 integrator_str: str = "Verlet", output_file: str = "data",
                 temperature_k: float = 293, friction: float = 0.01, potential_str: str = "LJ",
                 att_list: list = ["energy"],
                 pressure: float = 10e+6, compressibility: float = 10e-11):
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
        self.friction = float(friction)
        self.pressure = float(pressure)
        self.compressibility = compressibility
        self.potential = self.getPotential(potential_str)
        self.integrator = self.getIntegrator(integrator_str)
        self.attachments = self.getAttachment(att_list)

    @classmethod
    def initNVE(cls, temperature: float):
        return cls(temperature_k=temperature, integrator_str="NVE")

    @classmethod
    def initNVT(cls, temperature: float):
        return cls(temperature_k=temperature, integrator_str="NVT")

    @classmethod
    def initNPT(cls, temperature: float):
        return cls(temperature_k=temperature, integrator_str="NPT")

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
            return asap_EMT

        elif potential_lower in ["lj", "lennardjones", "lennard_jones"]:
            from ase.calculators.lj import LennardJones
            return LennardJones
        else:
            raise ValueError(f"Invalid potential function: {potential}")

    def getIntegrator(self, integrator: str):
        integrator_lower = integrator.lower()
        if integrator_lower in ["verlet", "nve"]:
            from ase.md.verlet import VelocityVerlet  # för NVE
            return functools.partial(VelocityVerlet, timestep=self.timestep)

        elif integrator_lower in ["langevin", "nvt"]:
            from ase.md.langevin import Langevin  # för NVT
            return functools.partial(Langevin, timestep=self.timestep, temperature_K=self.temperature_k,
                                     friction=self.friction)

        elif integrator_lower in ["berendsen", "npt"]:
            from ase.md.nptberendsen import NPTBerendsen
            return functools.partial(NPTBerendsen, timestep=self.timestep, temperature_K=self.temperature_k,
                                     pressure_au=self.pressure, compressibility_au=self.compressibility)

        else:
            raise ValueError(f"Invalid integrator: {integrator}")

    def getAttachment(self, attachments):
        pos_attachments = {'energy': self.printEnergy,
                           "momenta": self.printMomentum,
                           "center_of_mass": self.printCenterOfMass}
        for a in attachments:
            if a not in pos_attachments.keys():
                raise ValueError(f"Invalid attachment: {a}")
        return [pos_attachments[a] for a in attachments]

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

        dyn = self.integrator(atoms=atoms)

        material_name = str(atoms.symbols)
        print("MATERIALNAME: ", material_name)

        # save traj
        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms)

        for a in self.attachments:
            dyn.attach(functools.partial(a, atoms=atoms),
                       interval=self.interval)  # Attach the different functions for printing

        dyn.attach(traj.write, interval=self.interval)

        logger = MDLogger(dyn, atoms=atoms, logfile=f"{self.output_file}.log",
                          header=True, peratom=True, mode='a')  # Create a logger for writing data
        dyn.attach(logger, interval=self.interval)  # Attach logger

        dyn.run(self.steps)  # RUN

    def printEnergy(self, atoms):
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        etot = epot + ekin
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


if __name__ == "__main__":
    md = MDBase()
    md.runMD(atoms)
    vizualiseTraj("data.traj")
