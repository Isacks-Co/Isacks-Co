from copy import copy
class Frame:
    def __init__(self,data_dict):
        self.data_dict = data_dict
    
    def __getitem__(self, name):
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
