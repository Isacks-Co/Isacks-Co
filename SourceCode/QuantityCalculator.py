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

from Utils.unitConversions import auToGPascal, evToJ
from ase.eos import EquationOfState
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.units import kB


hbar = 0.6582119569

logger = logging.getLogger(__name__)


class QuantityCalculator:
    """
    Derived-quantity calculator for MD post-processing.

    All methods are implemented as ``@staticmethod`` functions and do not
    maintain internal state.

    Notes
    -----
    Callers are responsible for providing consistent units. This class does
    not validate units or sampling assumptions.
    """


    @staticmethod
    def computeSpecificHeatNVT(E_tot_seq, total_mass_amu, T):
        """
        Compute specific heat capacity from total energy fluctuations (NVT).

        Parameters
        ----------
        E_tot_seq : array-like
            Total energy time series (typically in eV).
        total_mass_amu : float
            Total mass of the system in atomic mass units (amu).
        T : float
            Temperature in Kelvin.

        Returns
        -------
        float
            Specific heat capacity in eV/(amu·K), assuming energies are in eV.

        Notes
        -----
        This uses a fluctuation expression based on the variance of the total
        energy time series.
        """

        E_tot_seq = np.array(E_tot_seq)

        e_mean = np.mean(E_tot_seq)
        e_2_mean = np.mean(np.array(E_tot_seq) ** 2)
        prefactor = 1 / (kB * T ** 2)
        specific_heat = prefactor * (e_2_mean - e_mean ** 2) / total_mass_amu  # Specific heat in ev/amu*K

        return specific_heat


    @staticmethod
    def computeSelfDiffusionCoefficient(msd_list, sample_spacing):
        """
        Estimate self-diffusion coefficient from mean-squared displacement (MSD).

        Uses the Einstein relation with an endpoint slope estimate:

            D ≈ (MSD_end - MSD_0) / (6 * (t_end - t_0))

        Parameters
        ----------
        msd_list : array-like
            Mean-squared displacement values (typically in Å^2).
        sample_spacing : float
            Time between MSD samples (typically in fs).

        Returns
        -------
        float
            Diffusion coefficient in Å^2/fs (if MSD is Å^2 and time is fs).

        Notes
        -----
        This uses only the first and last MSD points. For better statistics,
        a linear fit over a diffusive regime is often preferred.
        """
        # Find the actual elapsed time
        timestep_list = []
        for i in range(len(msd_list)):
            timestep = i * sample_spacing
            timestep_list.append([timestep, i])

        msd0 = msd_list[0]
        msd_final = msd_list[-1]
        t_0 = timestep_list[0][0]
        t_end = timestep_list[-1][0]
        D = (msd_final - msd0) / (t_end - t_0)

        return D / 6  # Å^2/fs

    @staticmethod
    def computeDebyeTemperature(V_A3, mass_u, N, G, K):
        """
        Estimate the Debye temperature from density and elastic moduli.

        Parameters
        ----------
        V_A3 : float
            Volume in Å^3.
        mass_u : float
            Total mass in atomic mass units (u).
        N : int
            Number of atoms.
        G : float
            Shear modulus (unit consistency required).
        K : float
            Bulk modulus (same unit convention as G).

        Returns
        -------
        float
            Debye temperature in Kelvin.

        Notes
        -----
        The factor ``/ 10.18`` is used as a unit conversion (comment suggests
        conversion related to fs/Å). Ensure `G`, `K`, and density are consistent
        with this convention.
        """
        rho = (mass_u / V_A3)

        transversal_sound_velocity = np.sqrt(G / rho)
        longitudinal_sound_velocity = np.sqrt((K + 4.0 * G / 3.0) / rho)
        sound_velocity = ((1.0 / 3.0) * (
                    1.0 / (longitudinal_sound_velocity ** 3) + 2.0 / (transversal_sound_velocity ** 3))) ** (-1.0 / 3.0)

        n = (N / V_A3)

        Theta_D = (hbar / kB) * ((6.0 * np.pi ** 2 * n) ** (
                    1.0 / 3.0)) * sound_velocity / 10.18

        return Theta_D

    @staticmethod
    def computeLindemannIndex(msd, nearest_neighbour):
        """
        Compute the Lindemann index.

        Parameters
        ----------
        msd : float
            Mean-squared displacement (typically in Å^2).
        nearest_neighbour : float
            Nearest-neighbour distance (typically in Å).

        Returns
        -------
        float
            Lindemann index (dimensionless).
        """
        return np.sqrt(msd) / nearest_neighbour

    @staticmethod
    def calculateModuli(C_matrix):
        """
        Compute isotropic elastic moduli from a stiffness matrix.

        Parameters
        ----------
        C_matrix : array-like of shape (6, 6)
            Stiffness matrix in Voigt notation.

        Returns
        -------
        tuple of float
            ``(bulk_modulus, shear_modulus, youngs_modulus)``.
        """
        bulk_modulus = (C_matrix[0, 0] + 2 * C_matrix[0, 1]) / 3
        G_shear = (C_matrix[3, 3] + C_matrix[4, 4] + C_matrix[5, 5] + C_matrix[1, 1] - C_matrix[0, 1]) / 5
        youngs_modulus = 9 * bulk_modulus * G_shear / (3 * bulk_modulus + G_shear)
        return bulk_modulus, G_shear, youngs_modulus

    @staticmethod
    def nearestNeighborsMean(atoms_sequence, start: int, end: int = None):
        """
        Compute mean nearest-neighbour distance over a range of configurations.

        For each configuration in the interval [start, end), builds an ASE
        neighbour list using ``natural_cutoffs`` and computes each atom's
        nearest-neighbour distance. The return value is the mean over all atoms
        and all configurations in the interval.

        Parameters
        ----------
        atoms_sequence : sequence of ase.Atoms
            Sequence of atomic configurations.
        start : int
            Start index (inclusive).
        end : int, optional
            End index (exclusive). If None, uses ``start + 1``.

        Returns
        -------
        float or None
            Mean nearest-neighbour distance in Å, or None if no neighbours
            could be found for some atom.
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
    def computeBulkModulus(stretch_sequence):
        """
        Fit an equation of state to compute the bulk modulus.

        Collects potential energy and volume from a stretch sequence, sorts by
        volume, and fits a Birch–Murnaghan EOS using ASE's
        :class:`ase.eos.EquationOfState`.

        Parameters
        ----------
        stretch_sequence : sequence
            Sequence of frames containing energy and volume information.

        Returns
        -------
        float
            Bulk modulus in GPa.

        Side Effects
        ------------
        Writes an EOS plot image file named ``Ag-eos.png``.

        Notes
        -----
        This function calls ``_get(frame, name)``, which must be provided by the
        surrounding codebase.
        """

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
