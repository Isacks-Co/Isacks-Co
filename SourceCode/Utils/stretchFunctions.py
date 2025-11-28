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
from ase.io.trajectory import Trajectory
from collections import defaultdict

logger = logging.getLogger(__name__)
from .unitConversions import auToGPascal, evToJ



def calculateCMatrix(strech_sequence):  # TODO Move so it is computed on the fly.

    betas = [[], [], [], [], [], []]
    for frame in strech_sequence:
        # Add the corresponding values to the right beta (strain direction)
        betas[frame.info["beta"]].append([frame.info["strain"], frame.info["stress"]])
    beta_arrays = [np.array(beta, dtype=object) for beta in betas]

    # Create dictionaries, one for each beta, and store the arrays of matrices with epsilon as key
    beta_dicts = [defaultdict(list) for i in range(6)]
    for i in range(6):
        for epsilon, matrix in beta_arrays[i]:
            beta_dicts[i][epsilon].append(matrix)

    averages = []
    for beta in beta_dicts:
        avg_data = []

        for epsilon, matrices in beta.items():
            # Take elementwise average over all the matrices for each epsilon
            stacked = np.stack(matrices)
            avg_matrix = stacked.mean(axis=0)
            avg_data.append(np.array((epsilon, avg_matrix), dtype=object))

        avg_data = sorted(avg_data, key=lambda x: x[0])
        averages.append(avg_data)
    C = np.zeros((6, 6))
    for i in range(6):
        # Line fit epsilon vs sigma to find each c_ij
        epsilons = np.array([x[0] for x in averages[i]], dtype=float)
        for j in range(6):
            sigmas = np.array([x[1][j] for x in averages[i]], dtype=float)
            C[j, i] = np.polyfit(epsilons, sigmas, 1)[0]
    C *= 160.21766208  # Convert to GPa
    return C


def _numericalC(stretch_sequence):  # TODO Move to do on the fly
    """
    Calculates the elastic constants C11, C22, C33, C12, C44

    Input
        ---

    Output
        C_from_U: matrix, where C11, C12, C44 = C[0,0], C[0,1], C[3,3]
    """
    betas = [[], [], [], [], [], []]
    # Prefer reconstructing 1D slices from 2D trajectory to avoid stale 1D data

    used_2d = False
    try:

        tol = 1e-12
        for frame in stretch_sequence:
            info = getattr(frame, 'info', {})
            b1 = info.get('beta1')
            b2 = info.get('beta2')
            if b1 is None or b2 is None:
                continue
            if int(b1) == int(b2) and np.isclose(float(info.get('strain1', 0.0)), 0.0,
                                                 atol=tol):  # Looks at difference between strains when creating 2D data
                i = int(b1)
                eps = float(info.get('strain2', 0.0))
                energy = float(info.get('total_energy', np.nan))
                stress = np.array(info.get('stress'), dtype=float)
                betas[i].append([eps, stress, energy])
        used_2d = any(len(b) > 0 for b in betas)
    except Exception as e2:
        logger.info(f"Failed")
    if not used_2d:
        # Fallback to legacy 1D if available
        try:
            stretch_trajectory = Trajectory(self.settings.output_file + "_stretch_data.traj")
            for frame in stretch_trajectory:
                energy = frame.info.get('total_energy')
                betas[int(frame.info['beta'])].append([float(frame.info['strain']), np.array(frame.info['stress'],
                                                                                             dtype=float),
                                                       float(energy)])
            logger.info("Used legacy 1D stretch trajectory to compute elastic constants.")
        except Exception as e:
            logger.info(f"No usable stretch data found (2D or 1D): {e}")
    beta_arrays = [np.array(beta, dtype=object) for beta in betas]
    beta_dicts = [defaultdict(list) for i in range(6)]
    averages = []

    for i in range(6):
        for eps, matrix, energy in beta_arrays[i]:
            beta_dicts[i][eps].append((matrix, energy))

    for beta in beta_dicts:
        avg_data = []

        for eps, matrices_and_energies in beta.items():
            # Take elementwise average for each epsilon
            matrices = np.stack([me[0] for me in matrices_and_energies])
            energies = np.array([me[1] for me in matrices_and_energies], dtype=float)
            avg_matrix = matrices.mean(axis=0)
            avg_energy = np.nanmean(energies)  # ignore possible NaNs
            avg_data.append(np.array((eps, avg_matrix, avg_energy), dtype=object))

        avg_data = sorted(avg_data, key=lambda x: x[0])
        averages.append(avg_data)

    second_deriv = np.zeros((6, 6))

    twoD_energies = stretch_sequence[-1].info["2D Energies"]
    strains_axis = stretch_sequence[-1].info["Strains axis"]
    number_of_pairs = stretch_sequence[-1].info["Number of pairs"]  # Should usually be 6

    for i in range(int(np.sqrt(number_of_pairs))):
        energy_1 = np.array([x[2] for x in averages[i]], dtype=float)
        for j in range(int(np.sqrt(number_of_pairs))):
            energy_2 = np.array([x[2] for x in averages[j]], dtype=float)
            if i == j:
                second_deriv[i, j] = secondOrderNumericalDerivative(strains_axis, [energy_1, energy_2])
                continue
            else:
                try:
                    second_deriv[i, j] = secondOrderNumericalDerivative(strains_axis, twoD_energies[i][j])
                except Exception as e:
                    logger.info(f"2D stretch calc failed for ({i},{j}), WHY THE FRICK???!: {e}")
                    second_deriv[i, j] = 0.0

    C_from_U = second_deriv / stretch_sequence[0].get_volume()
    logger.debug(f"C_from_U = \n {C_from_U * auToGPascal(1)} \n")
    logger.info(
        f" \n C_11 = {auToGPascal(C_from_U[0, 0])} \n C_12 = {auToGPascal(C_from_U[0, 1])} \n C_44 = {auToGPascal(C_from_U[3, 3])}")

    return C_from_U



def secondOrderNumericalDerivative(eps_list, U):
    """
    Approximate the mixed partial derivative ∂²U/(∂x ∂y) using available data.

    Supported inputs:
    - Case A: U is a 2D grid (n×n) of energies sampled on a common 1D grid `eps_list` along both axes.
      Then this function computes the central cross-difference mixed partial on interior points and returns its average.

    - Case B: U is a list/tuple of two 1D energy traces, U=[U1(x), U2(y)],
      each sampled on the same grid `eps_list`, with x == y.

    Returns:
        float: Average mixed partial estimate (or diagonal pure second in the identical-traces case).
    """
    x = list(eps_list).copy()
    y = list(eps_list).copy()

    # Try to detect a 2D grid first (most accurate if provided)
    try:
        import numpy as _np
        U_arr = _np.asarray(U)
        # logger.debug(f"U_arr = {U_arr}")
    except Exception:
        U_arr = None

    if U_arr is not None and getattr(U_arr, 'ndim', 0) == 2 and U_arr.shape[0] == U_arr.shape[1] == len(x):
        n = len(x)
        if n < 3:
            raise ValueError("At least 3 grid points are required for mixed partial central differences.")
        # Ensure numeric types
        x = _np.asarray(x, dtype=float)
        y = _np.asarray(y, dtype=float)
        U2D = _np.asarray(U_arr, dtype=float)
        d2 = _np.full_like(U2D, _np.nan, dtype=float)

        # Central difference approximation
        for i in range(1, n - 1):
            dx = x[i + 1] - x[i - 1]
            if dx == 0:
                raise ValueError(f"Non-unique x around index {i}.")
            for j in range(1, n - 1):
                dy = y[j + 1] - y[j - 1]
                if dy == 0:
                    raise ValueError(f"Non-unique y around index {j}.")
                d2[i, j] = (U2D[i + 1, j + 1] - U2D[i + 1, j - 1] - U2D[i - 1, j + 1] + U2D[i - 1, j - 1]) / (dx * dy)

        return float(_np.nanmean(d2[1:-1, 1:-1]))

    # Otherwise interpret U as two 1D traces
    if not isinstance(U, (list, tuple)) or len(U) != 2:
        raise ValueError("U must be either a 2D (n×n) grid or a list/tuple [U1, U2] of two 1D traces.")

    U1 = list(U[0]).copy()
    U2 = list(U[1]).copy()

    n = len(x)
    if not (len(y) == n and len(U1) == n and len(U2) == n):
        raise ValueError("x, y, U[0], and U[1] must all have the same length.")
    if n < 2:
        raise ValueError("At least two points are required to compute numerical derivatives.")

    # Convert to float
    try:
        x = [float(v) for v in x]
        y = [float(v) for v in y]
        U1 = [float(v) for v in U1]
        U2 = [float(v) for v in U2]
    except Exception as e:
        raise ValueError(f"All inputs must be numeric. Error: {e}")

    # Detect identical traces (diagonal case i == j in caller)
    if _np.allclose(U1, U2, rtol=1e-12, atol=1e-12):
        # Return mean pure second derivative along the shared axis
        _, d2 = numericalDerivative(x, U1, deg=2)
        return float(np.mean(d2))

    # logger.debug(f"returning 0.0")
    return 0.0