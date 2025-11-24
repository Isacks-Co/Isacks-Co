from copy import copy
import numpy as np

"""
Module used as a replacement of ASE Trajectories. 
This module only stores the desired data to sample
and not all the positions, forces .... that ASE did. 
Resulting in a much more memory efficient sampling. 
When reading data it works mostly the same as an ASE trajectory.
"""

class Frame:
    def __init__(self,time,data_dict):
        
        self.time = time
        self.data_dict = data_dict
    @property
    def keys(self):
        return self.data_dict.keys()
    
    @property
    def vals(self):
        return self.data_dict.values()
    def __getitem__(self, name):
        if name =="time":
            return self.time
        try:
            return self.data_dict[name]
        except KeyError:
            raise AttributeError(f"Frame has no attribute '{name}'")


class DataTrajectory:
    def __init__(self, initial_atomic_structure):
        self.initial_atoms = copy(initial_atomic_structure) 
        self._frames = [] 

    def append(self,frame: Frame):
        self._frames.append(frame)

    def storeTxtFile(self):
        col_width = 30 # Should work fine with current number of decimals
        with open(f"sampledata.txt", "w") as f:
            #HEADER
            f.write(f"{self.initial_atoms.label}\n")
            #f.write(f"Ensemble: {self.settings.ensemble}\n")
            
            #f.write(f"\n")

            #DATA
            
            f.write(f'{"time":<{col_width}}')
            f.write(f"".join(f"{label:<{col_width}}" for label in self._frames[0].keys) + "\n")
            for frame in self._frames[500:]:
                f.write(f"{frame.time:<{col_width}}")
                f.write(f"".join(f"{data:<{col_width}}" for data in frame.vals) + "\n")
                
            #f.write("".join(f"{value:<{col_width}.3f}" for value in quantities) + "\n")


    def __len__(self):
        return len(self._frames)
    
    def __iter__(self): 
        
        return iter(self._frames)
    
    def __getitem__(self, idx):
        if isinstance(idx, slice):
            new_traj = DataTrajectory(self.initial_atoms)
            new_traj._frames = self._frames[idx]
            return new_traj
        return self._frames[idx]
    
