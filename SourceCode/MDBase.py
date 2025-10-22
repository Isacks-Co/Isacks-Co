from simulationInput import SimulationSettings
import functools

import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary, ZeroRotation
from ase.units import fs, GPa
import logging
from LJRegistry import LJParams, calcMaxRc
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

        traj = Trajectory(filename=f"{self.output_file}.traj", mode="w", atoms=atoms) ## currently have .. before

        dyn.attach(lambda: self.save_data(atoms,traj),
                   interval=self.interval)

            # Add any other custom calculations here

        # Continue with the main MD run
        dyn.run(self.steps)  # RUN
        traj.close()  # Explicitly close the trajectory



        # Run stretch sequence for elastic constants
        if self.ensemble == "NVT":
            self._make_eos_traj(atoms)
            self._runStretchSequence(atoms)


    def save_data(self, atoms,traj):

        Ep = atoms.get_potential_energy()
        Ek = atoms.get_kinetic_energy()
        Et = atoms.get_total_energy()
        F = atoms.get_forces()
        V = atoms.get_volume()
        T = atoms.get_temperature()

        atoms.info['E_pot'] = float(Ep)
        atoms.info['E_kin'] = float(Ek)
        atoms.info['E_tot'] = float(Et)
        atoms.info['V'] = float(V)
        atoms.info['T'] = float(T)
        #atoms.arrays['F'] = F
        if 'F' in atoms.arrays:
            atoms.arrays['F'][:] = np.asarray(F, float)
        else:
            atoms.new_array('F', np.asarray(F, float))
        atoms.info['F'] = np.asarray(F, float).tolist()

        traj.write()


    def _runStretchSequence(self, atoms):

        def getStress(traj,atoms=atoms):
            atoms.info["stress"] = atoms.get_stress(voigt = True)
            atoms.get_total_energy()
            atoms.get_potential_energy()
            atoms.get_volume()
            traj.write()

        # Small strain amplitude
        eta = 5e-3  # 0.5%
        # Number of MD steps to run at each strained state
        hold_steps = 100
        I = np.eye(3)
        # Symmetric small-strain deformation gradients (F ≈ I + eps for small strains)
        F_list = []
        # ± isotropic
        F_list.append(I * (1.0 + eta))
        F_list.append(I * (1.0 - eta))
        # Orthorhombic (volume-conserving to first order): diag(1+eta, 1-eta, 1)
        #F_list.append(np.diag([1.0 + eta, 1.0 - eta, 1.0]))
        #F_list.append(np.diag([1.0 - eta, 1.0 + eta, 1.0]))
        # Symmetric shears: xy, xz, yz (F = I + eps, eps_ij = eps_ji = eta)
        eps_xy = I.copy()
        eps_xy[0, 1] = eps_xy[1, 0] = eta
        F_list.append(eps_xy)

        eps_xz = I.copy()
        eps_xz[0, 2] = eps_xz[2, 0] = eta
        F_list.append(eps_xz)

        eps_yz = I.copy()
        eps_yz[1, 2] = eps_yz[2, 1] = eta
        F_list.append(eps_yz)

        traj = Trajectory(filename=f"{self.output_file}_stretch_data.traj", mode="w", atoms=atoms)
        dyn = self.integrator(atoms=atoms)

        dyn.attach(lambda: getStress(traj=traj,atoms=atoms), 1)
        dyn.run(hold_steps)

        for F in F_list:
            A = atoms.cell.array.T
            A_new = (F @ A).T
            atoms.set_cell(A_new, scale_atoms=True)
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

    def _make_eos_traj(self, atoms, eta=0.12, npoints=20, relax=False, fmax=0.02,
                       traj_path=f"../Outputs/isotropic_stretch.traj"):
        """
        Creates isotrophic scaling of cell
            eta: max scaling factor
            npoints: nbr of points between (1-eta) och (1+eta)
            relax: if True, relax only inrternal positions at everyt volume
            fmax: force criteria for relax
            traj_path:
        Return: traj_path
        """
        from ase.optimize import BFGS
        from ase.md.langevin import Langevin


        if traj_path is None:
            traj_path = f"{self.output_file}_eos.traj"



        scales = np.linspace(1.0 - eta, 1.0 + eta, npoints)
        A0 = atoms.cell.array.copy()
        traj = Trajectory(traj_path, mode="w", )

        for s in scales:
            a = atoms.copy()
            a.set_cell(A0 * s, scale_atoms=True)  # isotrop scaling
            a.calc = self.potential(a)

            if hasattr(a.calc, "set_atoms"):
                a.calc.set_atoms(a)

            if relax:
                BFGS(a, logfile=None).run(fmax=fmax)  # relax internal coordi
            E_pot = a.get_potential_energy()
            V = a.get_volume()
            a.info['E_pot'] = float(E_pot)
            a.info['V'] = float(V)


            traj.write(a)
