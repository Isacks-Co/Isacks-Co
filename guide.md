# Project Overview

Small guide to run the program and set up the input files.

---

## Files in this project

- `MolecularDynamics.py` — Main file that runs the program.  
- `settings.json` —  Includes the parameters for the calculation
- `POSCAR` — Includes the atomic structure
 

---

## How to run the project

1. Open your terminal and navigate to the project folder. Make sure all the required files are in the same folder.  
2. Run the following command to start the calculation:

```bash
python MolecularDynamics.py -optionalflags
```

## Setting.json
The settings file consists of the following inputs:
```json
{ 
    "Temperature" : 300, 
    "Ensemble" : "NPT",
    "Potential" : "LJ",
    "Timestep" : 10,
    "Number_of_steps" : 1000,
    "Output_file" : "output", 
    "Friction" : 0,  
    "Compressibility" : 0,
    "Pressure" : 0,  
    "Interval" : 10, 
    "Supercells": 0 
}
```
### Parameters:
- Temperature -  As a positive float (K) * 
- Ensemble -  NVE, NVT, NPT *
- Potential -  LJ, EMT, MACE *
- Timestep - Time step as a float (fs) *
- Number_of_steps -  Total number of steps as an integer *
- Output_file - Name of the output file *
- Interval -  How often values are printed to the terminal *
- Supercells -  If you want a supercell the value 3 gives a 3 x 3 x 3 supercell. For a single unit cell use 0. * 
- Compressibility -  As a float
- Pressure -  As a float (Pa) 
- Friction  - As a float 


For all calculations the settings with (*) are required. For NPT, pressure and compressibility are also required. For NVT, friction is required. If a setting is missing.


## POSCAR
This file includes the atomic configuration. It follows the VASP standard. 
Below follows a minimal example.

```text
Cubic BN #Comment (e.g structure name)
3.57 # Lattice parameter (As a scaling factor on the lattice vectors)
0.0 0.5 0.5 # Lattice vectors
0.5 0.0 0.5
0.5 0.5 0.0
Cu N # Atomic species
1 1 # Number of each species
Direct # How the atomic positions are expressed (Direct - expressed in lattice vectors or Carthesian)
0.00 0.00 0.00 # Atomic positions  
0.25 0.25 0.25
```
The atomic positions is in order of the species, so in the example 0 0 0 is for the Copper atom and 0.25 0.25 0.25 is for the Nitrogen.
