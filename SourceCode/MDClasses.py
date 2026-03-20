# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import numpy as np
from ASEWrappers import AtomicStructure
from ASEWrappers import DataTrajectory, Frame
from Utils.equilibriumCondition import EquilibriumCondition
from ase.io.trajectory import Trajectory
from SimulationInput import SimulationSettings
from ase.units import fs
from copy import copy
import logging
log = logging.getLogger(__name__)
class MDBase:
    """
    Base class for molecular dynamics run types.

    Provides shared utilities for different MD workflows, including:
    - sampling scalar quantities from an AtomicStructure
    - storing sampled data in a DataTrajectory
    - attaching ASE trajectory writers to an integrator

    This class is not intended to be used directly, but to be subclassed
    by concrete MD run implementations.

    Attributes
    ----------
    integrator
        ASE-style integrator used to advance the simulation.
    sample_data : list or None
        Names of quantities to sample during the run.
    run_type : str or None
        Name of the run mode, used for output file naming.
    """
    def __init__(self, settings):
        self.integrator = settings.integrator
        self.sample_data = None
        self.run_type = None

    def _storeFrame(self, atomic_strucuture: AtomicStructure, data_traj: DataTrajectory, interval = 1):
        """
        Sample requested quantities and append a frame to a data trajectory.

        The simulation time is inferred from the previous frame time and the
        integrator timestep.

        Parameters
        ----------
        atomic_strucuture : AtomicStructure
            Current atomic configuration.
        data_traj : DataTrajectory
            Container storing sampled frames.
        """
        if len(data_traj) == 0:
            time = 0
        else:
            time = data_traj[-1].time + self.integrator.timestep * interval

        data = {label: self._getAtomsData(atomic_strucuture, label, data_traj.initial_atoms) for label in
                self.sample_data}

        frame = Frame(time, data)

        data_traj.append(frame)

    def _SaveASETrajectory(self, atomic_structure: AtomicStructure, interval=1):
        """
        Attach an ASE trajectory writer to the integrator.

        Parameters
        ----------
        atomic_structure : AtomicStructure
            Structure whose ASE Atoms object is written to disk.
        interval : int, optional
            Number of MD steps between trajectory frames.

        Notes
        -----
        The output file is named ``<run_type>.traj``.
        """
        traj = Trajectory(filename=f"{self.run_type}.traj", mode="w", atoms=atomic_structure.getAtoms())
        self.integrator.attach(traj.write, interval)

    @staticmethod
    def _getAtomsData(atomic_structure: AtomicStructure, name, initial_atomic_structure: AtomicStructure):
        """
        Extract a scalar quantity from an AtomicStructure.

        Parameters
        ----------
        atomic_structure : AtomicStructure
            Current structure state.
        name : str
            Name of the quantity to extract.
        initial_atomic_structure : AtomicStructure
            Reference structure used for displacement-based quantities.

        Returns
        -------
        float or object
            Extracted quantity corresponding to `name`.

        Raises
        ------
        RuntimeError
            If the structure has no associated potential.
        """

        if atomic_structure.potential is None:
            raise RuntimeError("Atoms object has no potential")

        if name == "E_pot":
            E_pot = atomic_structure.potential_energy
            return E_pot

        if name == "E_kin":
            E_kin = atomic_structure.kinetic_energy
            return E_kin

        if name == "E_tot":
            E_tot = atomic_structure.total_energy
            return E_tot

        if name == "V":
            vol = atomic_structure.volume
            return vol

        if name == "T":
            T = atomic_structure.temperature
            return T

        if name == "MSD":
            MSD = atomic_structure.computeMSD(initial_atomic_structure)
            return MSD

        if name == "P_internal":
            P_internal = atomic_structure.internal_pressure
            return P_internal

        if name == "NN":
            return atomic_structure.computeNearestNeighbour()


class EquilibriumRun(MDBase):
    """
    Equilibration molecular dynamics run.

    Runs the integrator until a fixed number of steps is reached or until
    an equilibrium condition is satisfied. The equilibrium criterion depends
    on the ensemble:
    - NVT: potential energy stability
    - other ensembles: internal pressure stability
    """
    def __init__(self, settings):
        super().__init__(settings)
        self.flag = 2
        self.run_type = "Equil"
        self.sample_data = ["T", "E_tot", "E_kin", "E_pot", "V", "MSD", "P_internal"]
        self.equil_data = np.array([])

    def run(self,atomic_structure: AtomicStructure ,num_steps,init_vel = True,store_traj = True, check_conv = False, check_expansion = False):
        """
        Run an equilibration simulation.

        Parameters
        ----------
        atomic_structure : AtomicStructure
            Structure to equilibrate (modified in-place).
        num_steps : int
            Maximum number of MD steps.
        init_vel : bool, optional
            If True, initialize velocities from a Maxwell–Boltzmann distribution.
        store_traj : bool, optional
            If True, write an ASE trajectory file.
        check_conv : bool, optional
            If True, check for equilibrium and stop early if satisfied.

        Returns
        -------
        AtomicStructure
            The equilibrated structure.
        """
        if init_vel:
            atomic_structure.setVelocitiesMB(self.integrator.temperature_K)

        if store_traj:
            self._SaveASETrajectory(atomic_structure, interval=500)

        if check_expansion:
            atoms_pre = atomic_structure.getAtoms()
            sorted_z_list = sorted([row[2] for row in atoms_pre.get_positions()])
            pre_height = sorted_z_list[-1] - sorted_z_list[0]
            self.integrator.attach(lambda: self._check_expansion_factor(atomic_structure, pre_height), 100)

        if check_conv:
            self.integrator.attach(lambda: self._saveData(atomic_structure),1)
            self.integrator.attach(lambda: self._check_equilibrium(atomic_structure),10)

        data_traj = DataTrajectory(atomic_structure)
        # Save first frame and sample every 100 interval
        self._storeFrame(atomic_structure, data_traj)
        self.integrator.attach(lambda: self._storeFrame(atomic_structure, data_traj, interval=500), 500)

        try:
            print("Starting equilibrium simulation")
            self.integrator.run(atomic_structure,num_steps)
        except StopIteration:
            print(f"Equilibrium reached in {len(self.equil_data)} steps")
        except Exception as err:
            print(err)

        data_traj.storeTxtFile(start_sample=0)
        with open("Output.txt", 'w') as o:
            o.write(f"{self.flag}\n{len(self.equil_data)}\n")
        return atomic_structure

    def _check_expansion_factor(self, atomic_structure : AtomicStructure, pre_height : float):
        """Stops simulation if the cell has doubled in height"""
        atoms_current = atomic_structure.getAtoms()
        sorted_z_list = sorted([row[2] for row in atoms_current.get_positions()])
        current_factor = sorted_z_list[-1] - sorted_z_list[0]
        expansion_factor = current_factor / pre_height
        if expansion_factor > 2:
            self.flag = 1
            print("Expansion factor exceeded 2, simulation stopped")
            raise StopIteration("Expansion factor exceeded 2, simulation stopped")

    def _check_equilibrium(self, atomic_structure):
        """
        Check whether equilibrium has been reached.

        Uses different stability criteria depending on the ensemble and
        raises ``StopIteration`` to terminate the MD run if equilibrium
        is detected.
        """
        if self.integrator.ensemble == "NVT":
            if len(self.equil_data) > 1000:
                if EquilibriumCondition.checkStable(self.equil_data[-1000:], 1e-6):
                    self.flag = 0
                    atomic_structure.final_energy_mean = np.mean(self.equil_data[-500:])
                    print(f"Stopped with equilibrium after {len(self.equil_data)}")
                    raise StopIteration(f"Equil reached")

        else:
            if len(self.equil_data) > 100:
                if EquilibriumCondition.checkInternalPressureStable(self.equil_data[-100:], 0.3):
                    self.flag = 0
                    print(f"Stopped with equilibrium after {len(self.equil_data)}")
                    raise StopIteration(f"Equil reached")


    def _saveData(self, atomic_structure):
        """
        Store equilibration diagnostic data.

        Depending on the ensemble, either potential energy or internal
        pressure is appended to the equilibration history.
        """
        if self.integrator.ensemble == "NVT":
            self.equil_data = np.append(self.equil_data, atomic_structure.total_energy)
        else:
            self.equil_data = np.append(self.equil_data, atomic_structure.internal_pressure * 160.21766208)
        
class SampleRun(MDBase):
    """
    Production MD run for sampling thermodynamic and structural quantities.

    Samples selected quantities at every MD step and stores them in a
    DataTrajectory, which is written to disk at the end of the run.
    """
    def __init__(self, settings : SimulationSettings):
        super().__init__(settings)
        self.run_type = "Sample"
        self.sample_data = ["T", "E_tot", "E_kin", "E_pot", "V", "MSD", "P_internal"]

        if settings.sample_nn:
            self.sample_data.append("NN")


    def run(self, atomic_structure: AtomicStructure, num_steps, store_traj=False):
        """
        Run a sampling simulation.

        Parameters
        ----------
        atomic_structure : AtomicStructure
            Structure to simulate (modified in-place).
        num_steps : int
            Number of MD steps to run.
        store_traj : bool, optional
            If True, write an ASE trajectory file.

        Returns
        -------
        None
        """
        data_traj = DataTrajectory(atomic_structure)
        if store_traj:
            self._SaveASETrajectory(atomic_structure)

        self.integrator.attach(lambda: self._storeFrame(atomic_strucuture=atomic_structure, data_traj=data_traj), 1)

        self.integrator.run(atomic_structure, num_steps)
        data_traj.storeTxtFile()
        return


class StretchRun(MDBase):
    """
    Small-strain MD run for estimating elastic stiffness coefficients.

    Applies a sequence of small strains, performs short MD holds for each
    strain, averages the resulting stress, and fits stress–strain slopes
    to construct a stiffness matrix.
    """
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Stretch"

    def run(self, atomic_structure: AtomicStructure):
        """
        Run a strain sweep and compute an elastic stiffness matrix.

        Parameters
        ----------
        atomic_structure : AtomicStructure
            Reference structure providing the initial cell and stress.

        Returns
        -------
        None

        Notes
        -----
        The resulting stiffness matrix is saved as ``cmatrix.npy``.
        """
        strains = np.linspace(-0.005, 0.005, 5)
        cell0 = atomic_structure.cell
        stress0 = atomic_structure.stress
        hold_steps = 25
        equil_atoms = copy(atomic_structure)

        C = np.zeros((6, 6))
        for beta in range(6):

            # list for storing the average matrix of stresses for each strain.
            average_stress = []
            for e in strains:
                # Strain tensor in Voigt form
                stress_list = []
                eps = np.zeros((3, 3))
                if beta < 3:
                    eps[beta, beta] = e
                elif beta == 3:
                    eps[1, 2] = eps[2, 1] = e / 2.0
                elif beta == 4:
                    eps[0, 2] = eps[2, 0] = e / 2.0
                elif beta == 5:
                    eps[0, 1] = eps[1, 0] = e / 2.0

                # Apply the strain to the cell and perform the number of steps specified with hold_steps
                new_cell = np.dot(cell0, np.eye(3) + eps)
                atoms = copy(equil_atoms)
                # atoms.calc = calculator
                atoms.cell = new_cell
                self.integrator.attach(lambda: self.appendStress(atoms, stress_list, stress0), 1)
                self.integrator.run(atoms, hold_steps)

                stacked = np.stack(stress_list)
                avg_matrix = stacked.mean(axis=0)
                average_stress.append(avg_matrix)
            sigmas = np.array(average_stress)
            for alpha in range(6):
                C[alpha, beta] = np.polyfit(strains, sigmas[:, alpha], 1)[0]
        np.save(f"cmatrix", C)
        return

    def appendStress(self, atoms, stress_list, stress0):
        """
        Append stress deviation relative to a reference stress.

        Parameters
        ----------
        atoms
            Structure providing a stress tensor.
        stress_list : list
            List to which stress values are appended.
        stress0
            Reference stress subtracted from the current stress.
        """
        # Help function for appending stresses during _stretchCell runs
        stress = atoms.stress - stress0

        stress_list.append(stress)
