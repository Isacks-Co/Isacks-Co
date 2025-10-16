from simulationInput import SimulationSettings
import functools
import logging

import numpy as np
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary, ZeroRotation
from ase.units import fs, GPa
import logging

from LJRegistry import LJParams, calcMaxRc # TODO Could these be squished into one function call?

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
            self.pressure = settings.pressure * GPa * 1e-9 # Pa to Au #TODO use functions
            self.compressibility = settings.compressibility/(GPa*1e-9) # Pa^-1 to Au
        
        #Integrator and potential
        self.integrator = self.getIntegrator(self.ensemble)
        self.potential = self.getPotential(settings.potential)



    def setupLJCalculator(self, atoms): # Can probably be moved. Maybe make a potential class??
        symbols = atoms.get_chemical_symbols()
        uniq = sorted(set(symbols))
        if len(uniq) != 1:
            raise ValueError(
                f"ASE LennardJones supports a single atom type only; found {uniq}. "
            )

        material_key = uniq[0].lower()
        params = LJParams(material=material_key)
        atomic_number = [(atoms.get_atomic_numbers()[0])] # TODO Add compatibility for multiple atom types
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
        equilibrium_steps = 10000
        dyn_eq = self.integrator(atoms=atoms)

         #dyn_eq.attach(lambda: self.checkConvergence(atoms), interval=max(1, int(1 / self.timestep)))
        #log.info(
          #  f"Starting equilibrium run with {self.ensemble} Ensemble to reach desired temperature of {self.temperature_k} K")
        #equil_traj = Trajectory(filename=f"../Outputs/equil_output_file.traj", mode="w", atoms=atoms) ## currently have .. before
        #dyn_eq.attach(lambda : self.compute_energies(atoms), interval = 20)
        #dyn_eq.attach(equil_traj.write, interval=20)
        try:
            dyn_eq.run(equilibrium_steps)

            current_T = atoms.get_temperature()
            log.info(f"Systems temperature is {round(current_T, 2)} K after {equilibrium_steps} steps")

        except RuntimeError as e:
            # pep 437 problem
            if "generator raised StopIteration" in str(e):
                log.info(
                    f"Equilibrium reached early (observer signaled StopIteration) at T = {round(atoms.get_temperature(), 2)}.") # ??
            else:
                raise RuntimeError(e)

        except StopIteration as ok:
            log.info(f"Equilibrium reached early: {ok}")

        except RuntimeWarning as err:
            log.warning(f"Equilibrium aborted due to instability: {err}")

        finally:
            self.quantity_list = []

    

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

        # Custom calculation function
        def save_custom_data(): # Should be moved
            """Store custom calculations in atoms.info"""
            atoms.info['stress eV/A3'] = atoms.get_stress(voigt=True)

            # Add any other custom calculations here

        dyn.attach(save_custom_data, interval=self.interval)

        #for a in self.attachments:
         #   dyn.attach(functools.partial(a, atoms=atoms),
         #              interval=self.interval)  # Attach the different functions for printing
        dyn.attach(lambda: self.save_data(atoms,traj),
                   interval=self.interval)
        #dyn.attach(traj.write, interval=self.interval)

        dyn.attach(lambda: self.checkDivergence(atoms),
                   interval=self.interval) 
        
        

        # Continue with the main MD run
        dyn.run(self.steps)  # RUN
        traj.close()  # Explicitly close the trajectory

    def save_data(self, atoms,traj):
        atoms.get_potential_energy()
        atoms.get_kinetic_energy()
        atoms.get_total_energy()
        atoms.get_forces()
        atoms.get_volume()
        atoms.get_positions()
        traj.write()

    def checkConvergence(self, atoms):
        self._updateQuantityList(atoms)
        if len(self.quantity_list) > 150 and not self.checkInstability(self.quantity_list):
            match self.ensemble:
                case "NVE":
                        raise StopIteration(f"Converged: E = {(atoms.get_potential_energy() + atoms.get_kinetic_energy()) / len(atoms):.2f} ev/atom")
                case "NVT":
                        raise StopIteration(f"Converged: T = {atoms.get_temperature():.2f} K")
                case "NPT":
                        raise StopIteration(f"Converged: T = {atoms.get_temperature():.2f} K")

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
        tol_energy_percentage = 0.05
        tol_temp = self.temperature_k * 0.10
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
                energy_tol = max(tol_energy_percentage * abs(mean_energy), 1e-4) # TODO Problematic since smaller enregy implies smaller margin for error
                nb_outside_tolerance = int(np.sum(np.abs(lastN - mean_energy) >= energy_tol))
                threshold = max(3, int(np.ceil(0.3 * num_points))) # 3 points must be outside of tolerance to trigger
                # Compares the number of points that are outside the threshold and returns True if they are

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
        #else:
        #    return (not nb_outside_tolerance > 10)
        #return False

    def _updateQuantityList(self, atoms):
        """
        Help function that appends quantities into quantity_list.
        """
        if not hasattr(self, "quantity_list"):
            self.quantity_list = [] #TODO Should not be class variable
        match self.ensemble:
            case "NVE":
                self.quantity_list.append((atoms.get_potential_energy() + atoms.get_kinetic_energy()) / len(atoms))
            case "NVT" | "NPT":
                self.quantity_list.append(atoms.get_temperature())
