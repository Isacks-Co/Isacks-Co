import numpy as np
import matplotlib.pyplot as plt
import os
import logging

logger = logging.getLogger(__name__)

def plot(x, y, save_as=None, resolution=None):
    # Fit once
    poly, error, x_fit, y_fit = fitter(x, y, deg=14)
    error_mean = np.mean(list(map(abs, error)))
    error_frac_list = error / error_mean

    """
    # High-res curve for top plot
    x_fit = np.linspace(x[0], x[-1], resolution) if resolution is not None else np.asarray(x)
    y_fit = poly(x_fit)

    # Residuals on original x
    x = np.asarray(x)
    y = np.asarray(y)
    y_pred_on_x = poly(x)
    error = np.abs(y - y_pred_on_x)  # or y - y_pred_on_x for signed residuals
    """

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

    # Top: data + fit
    ax1.axhline(0, color='gray', linewidth=1)
    ax1.scatter(x, y, s=7, label='Data', color='black')
    ax1.plot(x_fit, y_fit, label='Polyfit', color='blue', linestyle='--')
    ax1.set_xlabel("epsilon")
    ax1.set_ylabel("Data")
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