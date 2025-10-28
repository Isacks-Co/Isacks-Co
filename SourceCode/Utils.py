import numpy as np
import matplotlib.pyplot as plt

def plot(data, dt=1.0, save_as=None):
    
    data = np.array(data)
    time = np.arange(1, len(data)+1) * dt  # start from frame 1
    
    plt.figure(figsize=(6,4))
    plt.plot(time, data, marker='o', linestyle='-', color='tab:blue')
    plt.grid(True)
    plt.tight_layout()
    
    if save_as:
        plt.savefig(save_as, dpi=300)
        plt.close()
    else:
        plt.show()