import numpy as np
def plot(msd_list, dt=1.0, title="Mean Squared Displacement", save_as=None):
    import matplotlib.pyplot as plt

    """
    Plot MSD vs time.
    
    Parameters
    ----------
    msd_list : list or np.ndarray
        MSD values per frame.
    dt : float
        Timestep size in fs (default 1 fs).
    title : str
        Plot title.
    save_as : str or None
        Filename to save the plot. If None, plot is shown interactively.
    """
    msd_array = np.array(msd_list)
    time = np.arange(1, len(msd_array)+1) * dt  # start from frame 1
    
    plt.figure(figsize=(6,4))
    plt.plot(time, msd_array, marker='o', linestyle='-', color='tab:blue')
    plt.xlabel("Time [fs]")
    plt.ylabel("MSD [Å²]")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    
    if save_as:
        plt.savefig(save_as, dpi=300)
        plt.close()
    else:
        plt.show()