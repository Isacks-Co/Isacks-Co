import numpy as np
import functools
from ase.io.trajectory import Trajectory
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.units import fs, GPa
from ase.visualize import view
import logging
from simulationInput import SimulationSettings
from potentialSetUp import Potential

log = logging.getLogger(__name__)


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

        # General parameters
        self.ensemble = settings.ensemble
        self.timestep = settings.timestep * fs
        self.steps = settings.num_steps
        self.interval = settings.sample_interval
        self.output_file = settings.output_file

        # Ensemble specific parameters
        if settings.ensemble == "NVE":
            self.temperature_k = settings.initial_temperature

        elif settings.ensemble == "NVT":
            self.temperature_k = settings.temperature
            self.friction = settings.friction / fs

        elif settings.ensemble == "NPT":
            self.temperature_k = settings.temperature
            self.pressure = settings.pressure * GPa * 1e-9  # Pa to Au
            self.compressibility = settings.compressibility / (GPa * 1e-9)  # Pa^-1 to Au

        # Integrator and potential
        self.integrator = self.getIntegrator(self.ensemble)
        self.potential = Potential().getPotential(settings.potential)

    def getIntegrator(self, integrator: str):
        """
        Based on the argument return the integrator as a partial
        function with prefilled data
        args:
            str integrator: either NVE,NVT,NPT
        return:
            partial integrator: partial function of ASE integrator
        """
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



    def equilibriumRun(self, atoms):
        equilibrium_steps = 20000
        dyn_eq = self.integrator(atoms=atoms)

        dyn_eq.attach(lambda: self.checkConvergence(atoms), interval=max(1, int(1 / self.timestep)))
        # log.info(
        #  f"Starting equilibrium run with {self.ensemble} Ensemble to reach desired temperature of {self.temperature_k} K")

        try:
            dyn_eq.run(int(equilibrium_steps))
            current_T = atoms.get_temperature()
            log.info(f"Systems temperature is {round(current_T, 2)} K after {equilibrium_steps} steps")

        except RuntimeError as e:
            # pep 437 problem
            if "generator raised StopIteration" in str(e):
                log.warning(
                    f"Equilibrium reached early (observer signaled StopIteration) at T = {round(atoms.get_temperature(), 2)}.")
            else:
                raise RuntimeError(e)

        except StopIteration as ok:
            log.info(f"Equilibrium reached early: {ok}")

        except RuntimeWarning as err:
            log.error(f"Equilibrium aborted due to instability: {err}")
            quit()

        finally:
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

        # might be better to place these after eq run?
        Stationary(atoms, preserve_temperature=True)
        ZeroRotation(atoms, preserve_temperature=True)

        self.equilibriumRun(atoms=atoms)
        log.info("MD run starts with: %i steps", self.steps)
        dyn = self.integrator(atoms=atoms)

        traj = Trajectory(filename=f"../{self.output_file}.traj", mode="w", atoms=atoms)  ## currently have .. before

        # Custom calculation function
        def save_custom_data():
            """Store custom calculations in atoms.info"""
            atoms.info['stress eV/A3'] = atoms.get_stress(voigt=True)

            # Add any other custom calculations here

        dyn.attach(save_custom_data, interval=self.interval)

        # for a in self.attachments:
        #   dyn.attach(functools.partial(a, atoms=atoms),
        #              interval=self.interval)  # Attach the different functions for printing

        dyn.attach(traj.write, interval=self.interval)

        dyn.attach(lambda: self.checkDivergence(atoms),
                   interval=self.interval)

        # Apply a short sequence of slight, controlled strains and run a few steps at each.
        # This creates trajectory frames with non-zero strain for robust post-processing of elastic constants.
        def _apply_F_and_run(F, steps):
            A = atoms.cell.array.T
            A_new = (F @ A).T
            atoms.set_cell(A_new, scale_atoms=True)
            # log.info(f"Applied strain F=\n{F}\nCell now: {atoms.cell.cellpar()}")
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

        # log.info(f"Starting pre-production strain sequence with {len(F_list)} strains; {hold_steps} steps each")
        for F in F_list:
            _apply_F_and_run(F, hold_steps)

        # Continue with the main MD run
        dyn.run(self.steps)  # RUN
        traj.close()  # Explicitly close the trajectory

    def checkConvergence(self, atoms):
        """

        """
        self._updateQuantityList(atoms)
        window = 500
        error = 1e-6
        if len(self.quantity_list) < window:
            return

        last = np.asarray(self.quantity_list[-window:], dtype=float)
        mean_diff = abs(last[int(window/2):].mean() - last[:int(window/2)].mean())
        last_diff = abs(self.quantity_list[-1] - self.quantity_list[-2])
        match self.ensemble:
            case "NVT":
                error_per_atom = error*len(atoms)
                if last_diff <= error_per_atom and mean_diff <= error_per_atom:
                        raise StopIteration(f"Converged at E_pot = {last[-1]} eV.")

            case "NVE":
                pass

            case "NPT":
                raise StopIteration(f"Converged: T = {atoms.get_temperature():.2f} K")




    def _updateQuantityList(self, atoms):
        """
        Help function that appends quantities depending on the ensemble into quantity_list.
        """
        if not hasattr(self, "quantity_list"):
            self.quantity_list = []
        match self.ensemble:
            case "NPT":
                self.quantity_list.append(atoms.get_volume())
            case "NVT":
                self.quantity_list.append(atoms.get_potential_energy())
            case "NVE":
                self.quantity_list.append(atoms.get_temperature())
