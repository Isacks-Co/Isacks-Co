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
    Abstract class for MDRunner objects
    
    """
    def __init__(self, settings):
        self.integrator = settings.integrator
        self.sample_data = None
        self.run_type = None

    def _storeFrame(self, atomic_strucuture: AtomicStructure, data_traj: DataTrajectory):
        if len(data_traj) == 0:
            time = 0
        else:
            time = data_traj[-1].time + self.integrator.timestep

        data = {label: self._getAtomsData(atomic_strucuture, label, data_traj.initial_atoms) for label in
                self.sample_data}

        frame = Frame(time, data)

        data_traj.append(frame)

    def _SaveASETrajectory(self, atomic_structure: AtomicStructure, interval=1):
        """
        Function to save the ASE trajectory
        Args:
            atomic_structure (AtomicStructure): The atomic structure object containing our atoms, and other relevant
            information
            interval (int, optional): The interval between frames. Defaults to 1.
        """
        traj = Trajectory(filename=f"{self.run_type}.traj", mode="w", atoms=atomic_structure.getAtoms())
        self.integrator.attach(traj.write, interval)

    @staticmethod
    def _getAtomsData(atomic_structure: AtomicStructure, name, initial_atomic_structure: AtomicStructure):
        """
        Help function for getting specific data from the ASE atoms object.
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
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Equil"
        self.sample_data = ["E_pot"]
        self.equil_data = []

    def run(self,atomic_structure: AtomicStructure ,num_steps,init_vel = false,store_traj = True, check_conv = False):
        """
        The function that attaches other functions such as converge control etc and starts the EquilibriumRun simulation.
        Args:
            atomic_structure (AtomicStructure): The atomic stucture that we use in the simulation
            num_steps (int): The number of steps to run
            init_vel (bool, optional): Whether or not to initialize the velocity. Defaults to False.
            store_traj (bool, optional): Whether or not to store the trajectory. Defaults to True.
            check_conv (bool, optional): Whether or not to check the convergence. Defaults to False.
        """
        if init_vel:
            atomic_structure.setVelocitiesMB(self.integrator.temperature_K)

        if store_traj:
            self._SaveASETrajectory(atomic_structure)
        if check_conv:
            self.integrator.attach(lambda: self._saveData(atomic_structure),1)
            self.integrator.attach(self._check_equilibrium,10)
        try:
            self.integrator.run(atomic_structure,num_steps)
        except StopIteration as err:
            log.info(f"{err} in {len(self.equil_data)} steps")
        return atomic_structure
    
    def _check_equilibrium(self):
        """
        Function to check if the equilibrium condition is met.
        """
        if len(self.equil_data) > 100:
            if EquilibriumCondition.checkStable(self.equil_data[-100:], 0.01):
                print(f"Stopped with equilibrium after {len(self.equil_data)}")
                raise StopIteration(f"Equil reached")
    def _saveData(self, atomic_structure):
        """
        Function to save the data from the ASE atoms object.
        """
        self.equil_data.append(atomic_structure.potential_energy)

class SampleRun(MDBase):
    def __init__(self, settings : SimulationSettings):
        super().__init__(settings)
        self.run_type = "Sample"
        self.sample_data = ["T", "E_tot", "E_kin", "E_pot", "V", "MSD", "P_internal"]

        if settings.sample_nn:
            self.sample_data.append("NN")


    def run(self, atomic_structure: AtomicStructure, num_steps, store_traj=False):
        """
        The function that attaches other functions such as converge control etc and starts the SampleRun simulation.
        Args:
            atomic_structure (AtomicStructure): The atomic structure that we use in the simulation
            num_steps (int): The number of steps to run
            store_traj (bool, optional): Whether or not to store the trajectory. Defaults to False.

        """
        data_traj = DataTrajectory(atomic_structure)
        if store_traj:
            self._SaveASETrajectory(atomic_structure)

        self.integrator.attach(lambda: self._storeFrame(atomic_strucuture=atomic_structure, data_traj=data_traj), 1)

        self.integrator.run(atomic_structure, num_steps)
        data_traj.storeTxtFile()
        # TODO Add equil check / Fail check
        return


class StretchRun(MDBase):  # TODO Finish this
    def __init__(self, settings):
        super().__init__(settings)
        self.run_type = "Stretch"

    def run(self, atomic_structure: AtomicStructure):
        """
         The function that attaches other functions such as converge control etc and starts the Stretchrun simulation.
         Args:
            atomic_structure (AtomicStructure): The atomic structure that we use in the simulation

         """
        strains = np.linspace(-0.005, 0.005, 5)  # TODO Not hardcoded ?
        cell0 = atomic_structure.cell
        stress0 = atomic_structure.stress
        hold_steps = 500  # TODO Not hardcoded ?
        equil_atoms = copy(atomic_structure)
        calculator = atomic_structure.potential

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
        Function to append stress from the atoms to a list.
        Args:
            atoms (AtomicStructure): The atomic structure that we use in the simulation
            stress_list (list): The list containing stress values
            stress0 (float) : The initial stress
        """
        # Help function for appending stresses during _stretchCell runs
        stress = atoms.stress - stress0

        stress_list.append(stress)
