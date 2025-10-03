import functools
import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import fs,GPa
from ase.visualize import view
from SourceCode.logger import logger_setup
from SourceCode.simulationInput import SimulationSettings
from SourceCode.LJRegistry import LJParams, calcMaxRc




log = logger_setup()


class MDBase:
    """
    basic MD class
    When initialized it represents a MD-simulation with prefilled settings that can be used for multiple runs with
    different atomic structures.
    """
    def __init__(self, settings: SimulationSettings):
        """
        Create a MD runner from the settings specified in the 
        SimulationSettings object. 
        In:
            SimulationSettings object of type NVE,NVT or NPT
        
        """

        #General parameters
        self.ensemble = settings.ensemble
        self.timestep = settings.timestep * fs
        self.steps = settings.num_steps
        self.interval = settings.sample_interval
        self.output_file = settings.output_file
        

        #Ensemble specific parameters
        if settings.ensemble == "NVE":
            self.temperature_k = settings.initial_temperature

        elif settings.ensemble == "NVT":
            self.temperature_k = settings.temperature
            self.friction = settings.friction / fs

        elif settings.ensemble == "NPT":
            self.temperature_k = settings.temperature
            self.pressure = settings.pressure * GPa * 1e-9 # Pa to Au
            self.compressibility = settings.compressibility/(GPa*1e-9) # Pa^-1 to Au
        
        self.integrator = self.getIntegrator(self.ensemble)
        self.potential = self.getPotential(settings.potential)








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
        ro = params["ro_A"]
        rc = params["rc_A"]
        rc_max = calcMaxRc(atoms)

        if rc > rc_max:
            log.warning("The rCut is larger than the cell size, will use cell size to derive new value for rCut")
            rc = rc_max
            ro = 0.9*float(rc)

        print("using this rc ", rc)
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
                     material_key, atomic_number[0], eps, sig, rc)
            return calc_asap


        except Exception as e:
            from ase.calculators.lj import LennardJones as ase_LJ
            calc_ase = ase_LJ(
            epsilon=eps,
            sigma=sig,
            rc=rc,
            ro=ro
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
            from asap3.md.verlet import VelocityVerlet
            log.info("Integrator: Verlet")
            return functools.partial(VelocityVerlet, timestep=self.timestep)

        elif integrator_lower in ["langevin", "nvt"]:
            from asap3.md.langevin import Langevin
            log.info("Integrator: Langevin")
            return functools.partial(Langevin, timestep=self.timestep, temperature_K=self.temperature_k,
                                     friction=self.friction)

        elif integrator_lower in ["berendsen", "npt"]:
            from asap3.md.nptberendsen import NPTBerendsen
            log.info("Integrator: Berendsen")
            return functools.partial(NPTBerendsen, timestep=self.timestep, temperature_K=self.temperature_k,
                                     pressure_au=self.pressure, compressibility_au=self.compressibility)

        else:
            log.error("Invalid Integrator function: %s", integrator)  ##
            raise ValueError(f"Invalid integrator: {integrator}")

    def getAttachment(self, attachments):
        pos_attachments = {'energy': self.printEnergy,
                           "momenta": self.printMomentum,
                           "center_of_mass": self.printCenterOfMass,
                           "lattice": self.printLatticeConstants}

        for a in attachments:
            if a not in pos_attachments.keys():
                raise ValueError(f"Invalid attachment: {a}")


        return [pos_attachments[a] for a in attachments]

    def equilibriumRun(self, atoms):

        dyn_eq = self.integrator(atoms=atoms)
        self.equil_mode = True

        dyn_eq.attach(lambda: self.checkConvergence(atoms), interval=max(1, int(1 / self.timestep)))
        log.info(
            f"Starting equilibrium run with {self.ensemble} Ensemble to reach desired temperature of {self.temperature_k} K")

        try:
            dyn_eq.run(int(self.equilibrium_steps))

            current_T = atoms.get_temperature()
            log.info(f"Systems temperature is {round(current_T, 2)} K after {self.equilibrium_steps} steps")

        except RuntimeError as e:
            # Python 3.7+ translaterar StopIteration->RuntimeError inuti generatorer (PEP 479)
            if "generator raised StopIteration" in str(e):
                log.info(
                    f"Equilibrium reached early (observer signaled StopIteration) at T = {round(atoms.get_temperature(), 2)}.")
            else:
                raise RunTimeError(e)

        except StopIteration as ok:
            log.info(f"Equilibrium reached early: {ok}")

        except RuntimeWarning as err:
            log.warning(f"Equilibrium aborted due to instability: {err}")

        finally:
            self.equil_mode = False
            self.quantity_list = []

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

        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms)

        # Custom calculation function
        def save_custom_data():
            """Store custom calculations in atoms.info"""
            atoms.info['stress eV/A3'] = atoms.get_stress(voigt=True)

            # Add any other custom calculations here

        dyn.attach(save_custom_data, interval=self.interval)

        for a in self.attachments:
            dyn.attach(functools.partial(a, atoms=atoms),
                       interval=self.interval)  # Attach the different functions for printing

        dyn.attach(traj.write, interval=self.interval)

        logger = MDLogger(dyn, atoms=atoms, logfile=f"{self.output_file}.log",
                          header=True, peratom=True, mode='a')  # Create a logger for writing data
        dyn.attach(logger, interval=self.interval)  # Attach logger
        dyn.attach(lambda: self.checkDivergence(atoms),
                   interval=self.interval)  # TODO Possibly include checkConvergence here?

        # Apply a short sequence of slight, controlled strains and run a few steps at each.
        # This creates trajectory frames with non-zero strain for robust post-processing of elastic constants.
        def _apply_F_and_run(F, steps):
            A = atoms.cell.array.T
            A_new = (F @ A).T
            atoms.set_cell(A_new, scale_atoms=True)
            #log.info(f"Applied strain F=\n{F}\nCell now: {atoms.cell.cellpar()}")
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
        eps_xy = I.copy();
        eps_xy[0, 1] = eps_xy[1, 0] = eta;
        F_list.append(eps_xy)
        eps_xz = I.copy();
        eps_xz[0, 2] = eps_xz[2, 0] = eta;
        F_list.append(eps_xz)
        eps_yz = I.copy();
        eps_yz[1, 2] = eps_yz[2, 1] = eta;
        F_list.append(eps_yz)

        #log.info(f"Starting pre-production strain sequence with {len(F_list)} strains; {hold_steps} steps each")
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
        self._updateQuantityList(atoms)
        if len(self.quantity_list) > 150 and not self.checkInstability(self.quantity_list):
            if self.equil_mode:
                match self.ensemble:
                    case "NVE":
                            raise StopIteration(f"Converged: E = {(atoms.get_potential_energy() + atoms.get_kinetic_energy()) / len(atoms):.2f} ev/atom")
                    case "NVT":
                            raise StopIteration(f"Converged: T = {atoms.get_temperature():.2f} K")
                    case "NPT":
                            raise StopIteration(f"Converged: T = {atoms.get_temperature():.2f} K")

        # if len(quantity_list) > 1:
        #     recent = quantity_list[-window:]
        #     # Specifically for the equil run, to stop iterating when at the target temperature for ~5 iterations or more
        #     if self.ensemble != "NVE":
        #         lastN = recent[-7:]
        #         if (self.equil_mode and len(recent) >= 7):
        #             nb_in_tolerance = sum(abs(x - self.temperature_k) <= temperature_tol for x in lastN)
        #             if nb_in_tolerance >= 5:
        #                 raise StopIteration(
        #                     f"Target T reached: stopping at T = {atoms.get_temperature():.2f} K  "
        #                 )

    def checkDivergence(self, atoms):
        """
        Checks if temperature or energy diverges continously in one direction,
        returns an error if that's the case, uses a 'window' number of last runs to calculate mean average
        hits determines how many flucuations before exiting
        """
        # Uses helper function updateQuantityList.
        self._updateQuantityList(atoms)
        if (self.checkInstability(self.quantity_list)):
            # End the run and raise warning.
            raise RuntimeError(f"{self.ensemble}: Run cancelled.")

    def checkInstability(self, data_buffer):
        """
        Uses an input list to calculate stability, returns false if stable and true if unstable.

        Inputs: List with atomic temperature, energy or volume
        Output: Bool
        """
        window = 10
        if self.equil_mode:
            tol_temp = self.temperature_k * 0.05
        else:
            tol_temp = self.temperature_k * 0.10
        tol_energy_percentage = 0.05
        nb_outside_tolerance = 0
        threshold = 100000000 # needs to be large number
        if len(data_buffer) <= 1:
            return False
        match self.ensemble:
            case "NVT" | "NPT":
                lastN = data_buffer[-window:]
                nb_outside_tolerance = sum(abs(x - self.temperature_k) >= tol_temp for x in lastN)
            case "NVE":
                dt_eff_ps = self.interval * self.timestep / (1000 * fs)
                num_points = max(3, int(round(3.0 / dt_eff_ps))) # These two integers can be altered
                if len(data_buffer) < num_points:
                    return False
                lastN = np.asarray(data_buffer[-num_points:], dtype=float)
                mean_energy = float(np.mean(lastN))
                energy_tol = max(tol_energy_percentage * abs(mean_energy), 1e-4)
                nb_outside_tolerance = int(np.sum(np.abs(lastN - mean_energy) >= energy_tol))
                threshold = max(3, int(np.ceil(0.3 * num_points))) # 3 points must be outside of tolerance to trigger
                # Compares the number of points that are outside the threshold and returns True if they are
        if not self.equil_mode:
            if (nb_outside_tolerance >= threshold):
                log.warning(
                    f"Energy oscillates heavily — {nb_outside_tolerance}/{num_points} points outside tolerance"
                )
                return False
            if 12 < nb_outside_tolerance < 15:
                # Print warning in log, continues the run
                log.warning("Warning temperature oscillates.")
            if (nb_outside_tolerance >= 15):
                log.warning("Temperature oscillates heavily, cancelling run.")
                return True
        else:
            return (not nb_outside_tolerance > 10)
        return False

    def _updateQuantityList(self, atoms):
        """
        Help function that appends quantities into quantity_list.
        """
        if not hasattr(self, "quantity_list"):
            self.quantity_list = []
        match self.ensemble:
            case "NVE":
                self.quantity_list.append((atoms.get_potential_energy() + atoms.get_kinetic_energy()) / len(atoms))
            case "NVT" | "NPT":
                self.quantity_list.append(atoms.get_temperature())
