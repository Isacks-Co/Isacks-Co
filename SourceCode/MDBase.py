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
from SourceCode.LJRegistry import LJParams


EV_PER_A3_TO_GPA = 160.21766208
EV_PER_A3_TO_PA = 160.21766208e9  # Pa per (eV/Å^3)
AMU_TO_KG= 1.66053906892e-27
EV_TO_JOULE = 1.6021766208e-19
AVOGRADO = 6.02214076e23
A_TO_M = 1e-10

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

    def compressibilityAu(self, compressibility):
        compressibility_au = compressibility / 6.2415e-12
        return compressibility_au

    def setupLJCalculator(self, atoms):
        symbols = atoms.get_chemical_symbols()
        uniq = sorted(set(symbols))
        if len(uniq) != 1:
            raise ValueError(
                f"ASE LennardJones supports a single atom type only; found {uniq}. "
            )

        material_key = uniq[0].lower()
        params = LJParams(material=material_key)
        atomic_number = [(atoms.get_atomic_numbers()[0])]
        eps = params["epsilon_eV"]
        sig = params["sigma_A"]
        rc = params["rc_A"]


        try:
            from asap3 import LennardJones as asap_LJ
            calc_asap = asap_LJ(
                atomic_number,
                epsilon=[eps],
                sigma=[sig],
                rCut=rc,
                modified=True
            )

            atoms.calc = calc_asap
            _ = atoms.get_potential_energy()
            log.info("Using asap3 LJ | element=%s (Z=%s) | ε=%.4g eV | σ=%.4g Å | rc=%.4g Å ",
                     material_key, atomic_number[0], eps, sig, rc )
            return calc_asap


        except Exception as e:
            from ase.calculators.lj import LennardJones as ase_LJ
            calc_ase = ase_LJ(
            epsilon=eps,
            sigma=sig,
            rc=rc,
            ro=params["ro_A"]
            )
            log.warning(
                f"Falling back to ASE LJ | Reason: {e}"
            )
            return calc_ase






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
            return lambda atoms: asap_EMT()

        elif potential_lower in ["lj", "lennardjones", "lennard_jones"]:
            log.info("Potential: Lennard Jones")
            return self.setupLJCalculator
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

        gamma = self.friction * 3.0
        dyn_eq = self.integrator(atoms=atoms)
        real_ensemble  = getattr(self, "ensemble", None)
        self.equil_mode = True
        dyn_eq.attach(lambda: self.checkConvergence(atoms), interval= max(10, self.interval//2 ))
        log.info(f"Starting equilibrium run with {self.ensemble} Ensemble to reach desired temperature of {self.temperature_k} K")

        try:
            dyn_eq.run(int(self.equilibrium_steps))

            current_T = atoms.get_temperature()
            log.info(f"Systems temperature is {round(current_T,2)} K after {self.equilibrium_steps} steps")

        except RuntimeError as e:
            # Python 3.7+ translaterar StopIteration->RuntimeError inuti generatorer (PEP 479)
            if "generator raised StopIteration" in str(e):
                log.info(f"Equilibrium reached early (observer signaled StopIteration) at T = {round(atoms.get_temperature(), 2)}.")
            else:
                raise

        except StopIteration as ok:
            log.info(f"Equilibrium reached early: {ok}")

        except RuntimeWarning as err:
            log.warning(f"Equilibrium aborted due to instability: {err}")

        finally:
            if real_ensemble is not None:
               # self.ensemble = real_ensemble
                self.equil_mode = False
                self.temp_history = []




    def runMD(self, atoms):
        """
        In: 
            Atoms: ase Atoms object representing the crystal structure

        Runs a MD simulation with the setting specified in __init__
        Depending on attachments will possibly print some data.
        Will always save a trajectory and log file.        
        """
        atoms.calc = self.potential(atoms)
       
        MaxwellBoltzmannDistribution(atoms, temperature_K=self.temperature_k,
                                     force_temp=True)  # Initialize velocity according to temperature_k

        self.equilibriumRun(atoms=atoms)
        log.info("MD run starts with: %i steps", self.steps)
        dyn = self.integrator(atoms=atoms)

        #material_name = str(atoms.symbols)
        #print("MATERIALNAMN: ", material_name)

        # save traj
        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms)

        # Custom calculation function
        def save_custom_data():
            """Store custom calculations in atoms.info"""
            atoms.info['potential_energy eV'] = atoms.get_potential_energy()
            atoms.info['kinetic_energy eV'] = atoms.get_kinetic_energy()
            atoms.info['total_energy eV'] = atoms.get_total_energy()
            atoms.info['temperature'] = atoms.get_temperature()
            atoms.info['volume A3'] = atoms.get_volume()
            atoms.info['forces eV/A'] = atoms.get_forces()
            atoms.info['positions'] = atoms.get_positions()
            atoms.info['stress eV/A3'] = atoms.get_stress(voigt=True)
            atoms.info['number_of_atoms'] = atoms.get_global_number_of_atoms()
            atoms.info['cell'] = atoms.get_cell()
            atoms.info['cell_volume A3'] = atoms.get_cell().volume
            atoms.info['masses u'] = atoms.get_masses()
            atoms.info['density u/A3'] = sum(atoms.info['masses u']) / atoms.get_volume()


            # Add any other custom calculations here

        dyn.attach(save_custom_data, interval=self.interval)

        for a in self.attachments:
            dyn.attach(functools.partial(a, atoms=atoms),
                       interval=self.interval)  # Attach the different functions for printing

        dyn.attach(traj.write, interval=self.interval)

        logger = MDLogger(dyn, atoms=atoms, logfile=f"{self.output_file}.log",
                          header=True, peratom=True, mode='a')  # Create a logger for writing data
        dyn.attach(logger, interval=self.interval)  # Attach logger
        dyn.attach(lambda: self.failSafe(atoms), interval=self.interval)

        # Apply a short sequence of slight, controlled strains and run a few steps at each.
        # This creates trajectory frames with non-zero strain for robust post-processing of elastic constants.
        def _apply_F_and_run(F, steps):
            A = atoms.cell.array.T
            A_new = (F @ A).T
            atoms.set_cell(A_new, scale_atoms=True)
            log.info(f"Applied strain F=\n{F}\nCell now: {atoms.cell.cellpar()}")
            dyn.run(int(steps))

        # Small strain amplitude
        eta = 5e-3  # 0.5%
        # Number of MD steps to run at each strained state
        hold_steps = max(self.interval, 50)

        I = np.eye(3)
        # Symmetric small-strain deformation gradients (F ≈ I + eps for small strains)
        F_list = []
        # ± isotropic
        F_list.append(I * (1.0 + eta))
        F_list.append(I * (1.0 - eta))
        # Orthorhombic (volume-conserving to first order): diag(1+eta, 1-eta, 1)
        F_list.append(np.diag([1.0 + eta, 1.0 - eta, 1.0]))
        F_list.append(np.diag([1.0 - eta, 1.0 + eta, 1.0]))
        # Symmetric shears: xy, xz, yz (F = I + eps, eps_ij = eps_ji = eta)
        eps_xy = I.copy(); eps_xy[0,1] = eps_xy[1,0] = eta; F_list.append(eps_xy)
        eps_xz = I.copy(); eps_xz[0,2] = eps_xz[2,0] = eta; F_list.append(eps_xz)
        eps_yz = I.copy(); eps_yz[1,2] = eps_yz[2,1] = eta; F_list.append(eps_yz)

        log.info(f"Starting pre-production strain sequence with {len(F_list)} strains; {hold_steps} steps each")
        for F in F_list:
            _apply_F_and_run(F, hold_steps)

        # Continue with the main MD run
        dyn.run(self.steps)  # RUN
        traj.close()  # Explicitly close the trajectory

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


    def checkConvergence(self, atoms):
        window = 20
        temperature_tol = self.temperature_k * 0.05

        if (self.ensemble != "NVE"):
            self.temp_history.append(atoms.get_temperature())
        else:
            self.temp_history.append((atoms.get_potential_energy() + atoms.get_kinetic_energy()) / len(atoms))

        if len(self.temp_history) > 1:
            recent = self.temp_history[-window:]

            # Specifically for the equil run, to stop iteratign when at the target temperature for ~5 iterations or more
            if self.ensemble != "NVE":
                lastN = recent[-7:]
                if (self.equil_mode and len(recent) >= 7):
                    nb_in_tolerance = sum(abs(x - self.temperature_k) <= temperature_tol for x in lastN)
                    if nb_in_tolerance >= 5:
                        raise StopIteration(
                            f"Target T reached: stopping at T = {atoms.get_temperature():.2f} K  "
                        )


    def failSafe(self, atoms):
        """
        Checks if temperature diverges continously in one direction,
        returns an error if that's the case, uses a 'window' number of last runs to calculate mean average
        hits determines how many flucuations before exiting
        """
        window = 20
        tol_temp_percentage = 0.25
        tol_energy_percentage = 0.02

        if (self.ensemble != "NVE"):
            self.temp_history.append(atoms.get_temperature()) #temperature
        else:
            self.temp_history.append((atoms.get_potential_energy() + atoms.get_kinetic_energy()) / len(atoms)) #energy

        if len(self.temp_history) > 1:
            recent = self.temp_history[-window:]
            mean = np.mean(recent)
            std = np.std(recent)


            temperature_tol = self.temperature_k * tol_temp_percentage
            if self.ensemble != "NVE":
                if len(recent) >= 10:
                    lastN = recent[-10:]
                    nb_outside_tolerance = sum(abs(x - self.temperature_k) >= temperature_tol for x in lastN)
                    if nb_outside_tolerance >= 5:
                        #log.warning()
                        raise RuntimeWarning(
                            f"Run canceled because simulation is not stable. The temperature oscillations are greater than {tol_temp_percentage * 100}% of desired temperature.")

                    if self.ensemble == "NPT":
                        pass


            else:
                dt_ps = self.timestep / (1000.0 * fs)
                dt_eff_ps = self.interval * dt_ps
                window_ps = 2.0
                N = max(10, int(round(window_ps / dt_eff_ps)))

                if len(self.temp_history) >= N:
                    lastN = np.asarray(self.temp_history[-N:], dtype=float)
                    mean_energy = float(np.mean(lastN))
                    energy_tol = max(tol_energy_percentage * abs(mean_energy), 1e-4)
                    nb_outside = int(np.sum(np.abs(lastN - mean_energy) >= energy_tol))
                    thresh_hold = max(3, int(np.ceil(0.3 * N)))

                    if nb_outside >= thresh_hold:
                        raise RuntimeWarning(
                            f"NVE: unstable total energy/atom deviates more than "
                            f"{tol_energy_percentage * 100:.0f}% of the window mean "
                            f"(threshold = {energy_tol:.6e} eV/atom)"
                        )

                """
                 print(f"Standard devation: {std}, mean: {mean}")
                 if (abs(self.temperature_k - mean) > std)  or self.temp_history[-1] > 2*self.temperature_k:
                     # if std > self.temp_history[0]:
                     if self.hits == 5:
                         if (self.ensemble != "NVE"):

                         else:
                             raise RuntimeWarning("Run canceled because simulation is not stable. Total energy change is greater than 2 standard deviations.")
                     else:
                         self.hits += 1
                         print(self.hits)
                 """




                                