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
            self.tdamp = settings.tdamp
            self.pdamp = settings.pdamp


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

        elif integrator_lower in ["MKT", "npt"]:
            from asap3.md.nose_hoover_chain import IsotropicMTKNPT
            log.info("Integrator: Isotropic Martyna-Tobias-klein dynamics")
            return functools.partial(IsotropicMTKNPT, timestep=self.timestep, temperature_K=self.temperature_k,
                                     pressure_au=self.pressure, tdamp=self.tdamp, pdamp=self.pdamp)

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
        self._stretchCell(atoms)


    def save_data(self, atoms,traj):
        atoms.get_potential_energy()
        atoms.get_kinetic_energy()
        atoms.get_total_energy()
        atoms.get_forces()
        atoms.get_volume()
        atoms.get_positions()
        atoms.info["stress"] = atoms.get_stress()
        traj.write()

    def _runStretchSequence(self, atoms):

        # Small strain amplitude
        stretch_constant = 1e-2  # 1%
        # Number of MD steps to run at each strained state
        hold_steps = 1000
        stretch_steps = 5
        I = np.eye(3)
        # Symmetric small-strain deformation gradients (F ≈ I + eps for small strains)
        stretch_matrix_list = np.zeros((4*stretch_steps, 3, 3), dtype=float)
        type_list = np.empty(4*stretch_steps, dtype='<U10')
        count = 0

        for current_stretch in (np.linspace(-stretch_constant, stretch_constant, stretch_steps)):

            # Stretch in x direction
            stretch_xx = I
            stretch_xx[0,0] = 1 + current_stretch
            stretch_matrix_list[count] =stretch_xx
            type_list[count] = "stretch_xx"

            # Symmetric shears: xy, xz, yz (F = I + eps, eps_ij = eps_ji = eta)
            stretch_xy = I.copy()
            stretch_xy[0, 1] = stretch_xy[1, 0] = current_stretch
            stretch_matrix_list[stretch_steps + count] = stretch_xy
            type_list[stretch_steps + count] = "shears_xy"

            stretch_xz = I.copy()
            stretch_xz[0, 2] = stretch_xz[2, 0] = current_stretch
            stretch_matrix_list[2*stretch_steps + count] = stretch_xz
            type_list[2*stretch_steps + count] = "shears_xz"

            stretch_yz = I.copy()
            stretch_yz[1, 2] = stretch_yz[2, 1] = current_stretch
            stretch_matrix_list[3*stretch_steps + count] = stretch_yz
            type_list[3*stretch_steps + count] = "shears_yz"

            count += 1


        traj = Trajectory(filename=f"{self.output_file}_stretch_data.traj", mode="w", atoms=atoms)
        dyn = self.integrator(atoms=atoms)
        index  = 0
        dyn.attach(lambda: getStress(traj=traj,atoms=atoms, index=index), 1)

        def getStress(traj, atoms, index):
            atoms.info["stress"] = atoms.get_stress(voigt = True)
            atoms.info["stretch_matrix"]  = stretch_matrix_list[index]
            atoms.info["measurement"] = type_list[index]
            atoms.info["potential_energy"] = atoms.get_potential_energy()
            traj.write()

        dyn.run(1)

        for i in range(len(stretch_matrix_list)):
            index = i
            A = atoms.cell.array.T
            A_new = (stretch_matrix_list[i] @ A).T
            atoms.set_cell(A_new, scale_atoms=True)
            dyn.run(hold_steps)
            atoms.set_cell(A, scale_atoms=True)

    def elasticData(self, traj, atoms, strain, stress0, beta):
        atoms.info["strain"] = strain
        atoms.info["stress"] = atoms.get_stress(voigt=True) - stress0
        atoms.info["beta"] = beta
        traj.write()

    def _stretchCell(self, atoms):
        strains = np.linspace(-0.01, 0.01, 5)
        cell0 = atoms.get_cell()
        stress0 = atoms.get_stress(voigt=True)
        hold_steps = 500

        traj = Trajectory(filename=f"{self.output_file}_stretch_data.traj", mode="w", atoms=atoms)

        for beta in range(6):
            for e in strains:
                # Strain tensor in Voigt form
                eps = np.zeros((3, 3))
                if beta < 3:
                    eps[beta, beta] = e
                elif beta == 3:
                    eps[1, 2] = eps[2, 1] = e / 2.0
                elif beta == 4:
                    eps[0, 2] = eps[2, 0] = e / 2.0
                elif beta == 5:
                    eps[0, 1] = eps[1, 0] = e / 2.0

                # Apply the strain to the cell
                new_cell = np.dot(cell0, np.eye(3) + eps)
                atoms.set_cell(new_cell, scale_atoms=True)
                dyn = self.integrator(atoms=atoms)
                dyn.attach(lambda : self.elasticData(traj, atoms, e, stress0, beta))
                dyn.run(hold_steps)




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
