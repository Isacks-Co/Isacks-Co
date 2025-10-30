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


import copy
import logging
import matplotlib.pyplot as plt
import numpy as np
import os

logger = logging.getLogger(__name__)

def flattenDict(dict):

    flat = {}
    for key, value in dict.items():
        if isinstance(value, type(dict)):
            flat.update(flattenDict(value))
        else:
            flat[key] = value
    return flat



def plot(x, y, x_label="x", y_label="y", save_as=None):
    # Fit once
    poly, error, x_fit, y_fit = fitter(x, y, deg=2)
    error_mean = np.mean(list(map(abs, error)))
    error_frac_list = error / error_mean

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

    # Top: data + fit
    ax1.axhline(0, color='gray', linewidth=1)
    ax1.scatter(x, y, s=7, label=y_label, color='black')
    ax1.plot(x_fit, y_fit, label='Polyfit', color='blue', linestyle='--')
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax1.legend()

    # Bottom: absolute error vs x
    ax2.axhline(0, color='gray', linewidth=1)
    ax2.plot(x, error_frac_list, label='Absolute error ratio |y - y_fit(x)|/mean_error', color='red')
    ax2.set_xlabel("epsilon")
    ax2.set_ylabel("Error")
    ax2.legend()

    fig.suptitle(os.path.basename(save_as) if save_as else "Plot")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save_as:
        folder = os.path.dirname(save_as) or "../Outputs/Plots"
        os.makedirs(folder, exist_ok=True)
        full_path = os.path.join(folder, os.path.basename(save_as))
        fig.savefig(full_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()


def fitter(x, y, deg=14, resolution=None):
    coeffs = np.polyfit(x, y, deg=deg)
    poly = np.poly1d(coeffs)

    x_fit = np.linspace(x[0], x[-1], resolution) if resolution is not None else np.asarray(x)
    y_fit = poly(x_fit)

    # Residuals on original x
    x = np.asarray(x)
    y = np.asarray(y)
    y_pred_on_x = poly(x)
    error = np.abs(y - y_pred_on_x)  # or y - y_pred_on_x for signed residuals

    return poly, error, x_fit, y_fit


def numericalDerivative(x, y, deg=1):
    """
    Compute the numerical derivative dy/dx given lists of x and y values.
    Uses central difference for interior points and forward/backward for endpoints.

    Parameters:
        x (list or array): x-values (must be increasing and same length as y)
        y (list or array): y-values

    Returns:
        list: derivative values (same length as x)
    """
    iteration = 0
    while iteration < deg:
        if len(x) != len(y):
            raise ValueError("x and y must have the same length.")
        if len(x) < 2:
            raise ValueError("At least two points are required to compute a derivative.")
        dydx = [0.0] * len(x)
        x_out = [0.0] * len(x)

        # Forward difference for first point
        dydx[0] = (y[1] - y[0]) / (x[1] - x[0])
        x_out[0] = (x[1] + x[0]) / 2

        # Central difference for interior points
        for i in range(1, len(x) - 1):
            if x[i + 1] != x[i - 1]:
                dx = x[i + 1] - x[i - 1]
                dy = y[i + 1] - y[i - 1]
                # logger.debug(f"\n [{i}] ----- dx = {dx} \n [{i}] ----- dy = {dy} \n")
                dydx[i] = dy / dx
                x_out[i] = x[i]
            else:
                dydx[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                x_out[i] = (x[i + 1] + x[i]) / 2
                logger.info(f"Warning: x values are not unique at index {i}, using central difference.")
        # Backward difference for last point
        dydx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
        x_out[-1] = (x[-1] + x[-2]) / 2
        y = dydx
        x = x_out
        iteration += 1

    return x_out, dydx


def sortRealignAndFilter(x_list, y_list, filter=True, resolution=None):
    x_sorted = sorted(x_list)
    y_realigned = []
    for x in x_sorted:
        y_realigned.append(y_list[x_list.index(x)])

    poly, error, x_fit, y_fit = fitter(x_sorted, y_realigned, deg=2, resolution=resolution)
    error_mean = np.mean(list(map(abs, error)))

    if filter:
        to_remove = []
        for j, (quantity, index) in enumerate(zip(y_realigned, x_sorted)):
            if (error[j] / error_mean) >= 2.0:
                logger.debug(f"error ratio at step {j} : {error[j] / error_mean})")
                to_remove.append(j)

        for i in sorted(to_remove, reverse=True):
            print(f"Removing index {i} -> {y_realigned[i]}")
            del y_realigned[i]
            del x_sorted[i]
    return x_sorted, y_realigned


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
