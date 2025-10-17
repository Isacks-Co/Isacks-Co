from simulationInput import SimulationSettings
import functools

import numpy as np
import functools
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary, ZeroRotation
from ase.units import fs, GPa
import logging
from simulationInput import SimulationSettings
from potentialSetUp import Potential

log = logging.getLogger(__name__)


class MDBase: #TODO Look at unit conversions
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



    def equilibriumRun(self, atoms, equilibrium_steps = 40000):
        """
        Function that runs certain amount of steps 'equilibrium_steps', to make the system reach equilibrium,
        before the 'real' production simulation is run.

        in:
        atoms : ase.Atoms object
        """

        dyn_eq = self.integrator(atoms=atoms)
        dyn_eq.attach(lambda: self._checkDivergence(atoms), interval=max(1, int(1 / self.timestep) ))

        try:
            dyn_eq.run(equilibrium_steps)
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

        #finally:
        #    self.quantity_list = []

    

    def runMD(self, atoms): #TODO Needs better comments
        """
        In: 
            Atoms: ase Atoms object representing the crystal structure

        Runs a MD simulation with the setting specified in __init__
        Depending on attachments will possibly print some data.
        Will always save a trajectory and log file.        
        """
        atoms.calc = self.potential(atoms) # Still dont like this

        MaxwellBoltzmannDistribution(atoms, temperature_K=self.temperature_k,
                                     force_temp=True)  # Initialize velocity according to temperature_k
        Stationary(atoms) # Make sure center of mass has no linear momentum
        ZeroRotation(atoms) # Make sure center of mass has no angular momentum, might not be needed
        self.equilibriumRun(atoms=atoms) # TODO BREAKS TO EARLY

        log.info("MD run starts with: %i steps", self.steps)
        dyn = self.integrator(atoms=atoms)
        #dyn.attach(lambda: self._checkDivergence(atoms), interval=max(1, int(1 / self.timestep) ))

        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms) ## currently have .. before

        dyn.attach(lambda: self.save_data(atoms,traj),
                   interval=self.interval)

        #dyn.attach(lambda: self._checkDivergence(atoms),
                 #  interval=self.interval)

        # Continue with the main MD run
        dyn.run(self.steps)  # RUN
        traj.close()  # Explicitly close the trajectory

        # Run stretch sequence for elastic constants
        self._runStretchSequence(atoms)


    def save_data(self, atoms,traj):
        atoms.get_potential_energy()
        atoms.get_kinetic_energy()
        atoms.get_total_energy()
        atoms.get_forces()
        atoms.get_volume()
        atoms.get_positions()
        traj.write()

    def _runStretchSequence(self, atoms):

        # Small strain amplitude
        stretch_constant = 1e-2  # 1%
        # Number of MD steps to run at each strained state
        hold_steps = 100
        I = np.eye(3)
        # Symmetric small-strain deformation gradients (F ≈ I + eps for small strains)
        stretch_matrix_list = []
        measurement_list = []
        # Reference measurement
        stretch_matrix_list.append(I)
        measurement_list.append("reference")

        # ± isotropic
        stretch_matrix_list.append(I * (1.0 + stretch_constant))
        stretch_matrix_list.append(I * (1.0 - stretch_constant))
        measurement_list.append("isotropic_plus")
        measurement_list.append("isotropic_minus")

        # Orthorhombic (volume-conserving to first order): diag(1+eta, 1-eta, 1)
        stretch_matrix_list.append(np.diag([1.0 + stretch_constant, 1.0 - stretch_constant, 1.0]))
        stretch_matrix_list.append(np.diag([1.0 - stretch_constant, 1.0 + stretch_constant, 1.0]))
        measurement_list.append("orthorhombic_plus_minus")
        measurement_list.append("orthorhombic_minus_plus")

        # Symmetric shears: xy, xz, yz (F = I + eps, eps_ij = eps_ji = eta)
        stretch_xy = I.copy()
        stretch_xy[0, 1] = stretch_xy[1, 0] = stretch_constant
        stretch_matrix_list.append(stretch_xy)
        measurement_list.append("shears_xy")

        stretch_xz = I.copy()
        stretch_xz[0, 2] = stretch_xz[2, 0] = stretch_constant
        stretch_matrix_list.append(stretch_xz)
        measurement_list.append("shears_xz")

        stretch_yz = I.copy()
        stretch_yz[1, 2] = stretch_yz[2, 1] = stretch_constant
        stretch_matrix_list.append(stretch_yz)
        measurement_list.append("shears_yz")

        traj = Trajectory(filename=f"{self.output_file}_stretch_data.traj", mode="w", atoms=atoms)
        dyn = self.integrator(atoms=atoms)

        dyn.run(hold_steps)

        for i in range(len(stretch_matrix_list)):
            def getStress(traj,atoms=atoms, index=i):
                atoms.info["stress"] = atoms.get_stress(voigt = True)
                atoms.info["stretch_matrix"]  = stretch_matrix_list[index]
                atoms.info["measurement"] = measurement_list[index]
                atoms.info["steps"] = hold_steps
                traj.write()

            A = atoms.cell.array.T
            A_new = (stretch_matrix_list[i] @ A).T
            atoms.set_cell(A_new, scale_atoms=True)
            dyn.attach(lambda: getStress(traj=traj,atoms=atoms, index=i), 1)
            dyn.run(hold_steps)
            atoms.set_cell(A, scale_atoms=True)


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

    def _checkDivergence(self, atoms):
        """
        Very simple divergence check, that controls that the temperature is finite and that it does not
        rise above a threshold based on the variable 'div_factor. Raises errors if we pass it or if temperature
        becomes infinite.

        In:
            Atoms: ase Atoms object representing the crystal structure

        """

        div_factor = 1.8
        current_T = atoms.get_temperature()

        if not np.isfinite(current_T):
            raise RuntimeError(f"Divergence appeared, temperature is NaN/Inf")

        if (current_T > div_factor * self.temperature_k and self.temperature_k != 0):
            raise RuntimeError(f"Divergence: T={current_T:.1f} K >  {div_factor * self.temperature_k:.1f} K, "
                               f"(div_factor * desired temp).")


    def _updateQuantityList(self, atoms):
        """
        Help function that appends quantities depending on the ensemble into quantity_list.
        """
        if not hasattr(self, "quantity_list"):
            self.quantity_list = [] #TODO Should not be class variable
        match self.ensemble:
            case "NPT":
                self.quantity_list.append(atoms.get_volume())
            case "NVT":
                self.quantity_list.append(atoms.get_potential_energy())
            case "NVE":
                self.quantity_list.append(atoms.get_temperature())
