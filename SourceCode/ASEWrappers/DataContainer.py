from copy import copy
import numpy as np
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
            
            f.write(f"{"time":<{col_width}}")
            f.write(f"".join(f"{label:<{col_width}}" for label in self._frames[0].keys) + "\n")
            for frame in self._frames[500:]:
                f.write(f"{frame.time:<{col_width}}")
                f.write(f"".join(f"{data:<{col_width}}" for data in frame.vals) + "\n")
                
            #f.write("".join(f"{value:<{col_width}.3f}" for value in quantities) + "\n")


    @classmethod
    def fromTxtFile(cls, filename, initial_atomic_structure):
        """Rebuild a DataTrajectory from a text file."""
        with open(filename, "r") as f:
            lines = f.readlines()

   
        headers = lines[0].split()


        traj = cls(initial_atomic_structure)

        for line in lines[1:]:
            if not line.strip():
                continue 
            values = line.split()
         
            data_dict = {}
            for key, val in zip(headers, values):
                try:
                    data_dict[key] = float(val)
                except ValueError:
                    data_dict[key] = val
            traj.append(Frame(data_dict))

        return traj

    def __len__(self):
        return len(self._frames)
    
    def __iter__(self): 
        
        return iter(self._frames)
    
    def __getitem__(self, idx):
        if isinstance(idx, slice):
            new_traj = DataTrajectory()
            new_traj._frames = self._frames[idx]
            return new_traj
        return self._frames[idx]
