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


import logging
import numpy as np
from Utils.plotting import secondOrderNumericalDerivative
from Utils.unitConversions import auToGPascal, evToJ
from ase.eos import EquationOfState
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.units import kB
from scipy.constants import physical_constants

hbar = physical_constants['Planck constant over 2 pi in eV s'][0] * 1e15

logger = logging.getLogger(__name__)


class QuantityCalculator:
    """
    Computes derived quantities from sequences of instantaneus data
    """

    @staticmethod
    def computeSpecificHeatNVT(E_tot_seq, total_mass_amu, T):
        """ 
        Compute specific heat capacity as a time average of the total energy fluctuations 
        over all frames.
        Unit: ev/(amu*K) (amu = atomic mass unit)
        """

        E_tot_seq = np.array(E_tot_seq)

        e_mean = np.mean(E_tot_seq)
        e_2_mean = np.mean(np.array(E_tot_seq) ** 2)
        prefactor = 1 / (kB * T ** 2)
        specific_heat = prefactor * (e_2_mean - e_mean ** 2) / total_mass_amu  # Specific heat in ev/amu*K

        return specific_heat

    @staticmethod
    def computeSpecificHeatNVE(E_kin_seq, total_mass_amu, T):
        """
        Compute specific heat capacity as a time average of the kinetic energy fluctuations 
        over all frames. 
        Unit: ev/(amu*K) (amu = atomic mass units)
        """

        e_kin_mean = np.mean(E_kin_seq)
        e_kin_2_mean = np.mean(np.array(E_kin_seq) ** 2)

        specific_heat = (3 * kB / 2) * 1 / (
                    1 - (2 / (3 * (kB * T) ** 2) * (e_kin_2_mean - e_kin_mean ** 2))) / total_mass_amu

        return specific_heat

    @staticmethod
    def computeSelfDiffusionCoefficient(msd_list,
                                        sample_spacing):  # Needs constant temperature, for current implementation.

        """
        Compute self diffusion coefficient from the slope of MSD over a large timeperiod
        Unit: Å^2/fs
        """

        # Find the actual elapsed time
        timestep_list = []
        for i in range(len(msd_list)):
            timestep = i * sample_spacing
            timestep_list.append([timestep, i])

        if len(timestep_list) > 100:  # Since early values of MSD are inaccurate
            msd0 = msd_list[50]
            msd_final = msd_list[-1]
            t_0 = timestep_list[50][0]
            t_end = timestep_list[-1][0]
            D = (msd_final - msd0) / (t_end - t_0)

        else:
            logger.error("Too small sample size to calculate self-diffusion coefficient")
            D = None

        return D / 6  # Å^2/fs

    @staticmethod
    def computeDebyeTemperature(V_A3, mass_u, N, G, K):

        rho = (mass_u / V_A3)

        transversal_sound_velocity = np.sqrt(G / rho)
        longitudinal_sound_velocity = np.sqrt((K + 4.0 * G / 3.0) / rho)
        sound_velocity = ((1.0 / 3.0) * (
                    1.0 / (longitudinal_sound_velocity ** 3) + 2.0 / (transversal_sound_velocity ** 3))) ** (-1.0 / 3.0)

        n = (N / V_A3)

        Theta_D = (hbar / kB) * ((6.0 * np.pi ** 2 * n) ** (
                    1.0 / 3.0)) * sound_velocity / 10.18  # TODO move to unit conversion file  to fs/Å

        return Theta_D

    @staticmethod
    def computeLindemannIndex(msd, nearest_neighour_d):

        return np.sqrt(msd) / nearest_neighour_d

    @staticmethod
    def calculateModuli(C_matrix):
        bulk_modulus = (C_matrix[0, 0] + 2 * C_matrix[0, 1]) / 3
        G_shear = (C_matrix[3, 3] + C_matrix[4, 4] + C_matrix[5, 5] + C_matrix[1, 1] - C_matrix[0, 1]) / 5
        youngs_modulus = 9 * bulk_modulus * G_shear / (3 * bulk_modulus + G_shear)
        return bulk_modulus, G_shear, youngs_modulus

    # TODO Still need to look at the stuf below this

    @staticmethod
    def nearestNeighborsMean(atoms_sequence, start: int,
                             end: int = None):  # TODO Look at this and make work without traj
        """Calculate the mean distance of nearest neighbor in the structure for the last ten states of the simulation
        Loop structure: Last ten states -> Each atom -> neighbors to current atom

        (int) start : Index for the start of the interval that should be checked
        (int) end : Index for the end of the interval that should be checked
        """
        if end is None:
            end = start + 1
        INF = 1e9
        NN_list = []

        for state in range(start, end):
            # Load neighbor list for the current state
            atoms = atoms_sequence[state]
            cutoff = natural_cutoffs(atoms)
            neighbor_list = NeighborList(cutoff, bothways=True)
            neighbor_list.update(atoms)

            for current_atom in range(atoms.get_global_number_of_atoms()):
                # Loop over all atoms in and find their nearest neighbor
                indices, offsets = neighbor_list.get_neighbors(current_atom)
                nearest_distance = INF

                # First object seems to be the atom itself, don't loop over it
                for neighbor_index, offset in zip(indices[1:], offsets[1:]):
                    # Create a vector between current_atom and the neighbors in the list, save the shortest distance
                    NN_vector = atoms.positions[neighbor_index] + offset @ atoms.get_cell() - atoms.positions[
                        current_atom]
                    distance = np.sqrt(NN_vector.dot(NN_vector))

                    if distance < nearest_distance:
                        nearest_distance = distance

                if nearest_distance == INF:
                    error_msg = f"Could not calculate NN distance, didn't find any NN for atom {current_atom}"
                    logger.error(error_msg)
                    return

                NN_list.append(nearest_distance)
        NN_mean_distance = np.mean(NN_list)
        logger.debug(f"Mean value of nearest neighbor : {NN_mean_distance} å")
        return NN_mean_distance

    @staticmethod
    def computeBulkModulus(stretch_sequence):  # TODO Move so it is computed on the fly.
        energies = []
        cells = []
        for frame in stretch_sequence:
            energies.append(_get(frame, "E_pot"))
            cells.append(_get(frame, "V"))

        V = np.array(cells)
        E = np.array(energies)
        order = np.argsort(V)
        E, V = E[order], V[order]

        eos = EquationOfState(V, E, eos='birchmurnaghan')
        v0, e0, B0_eVa3 = eos.fit()
        B0_GPa = auToGPascal(B0_eVa3)
        eos.plot('Ag-eos.png')
        return B0_GPa
