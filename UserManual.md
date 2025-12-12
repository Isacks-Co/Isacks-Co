# User guide 📠
Below is the guide on how to setup and run Isack-Co's Molecular dynamic program.

## External dependencies 🎰
To run the program a set external modules are needed. To get these 
1. Navigate to the Documents folder.
2. run 
```bash
      pip install requirements.txt
```

## How to run the project  🚀
1. Open your terminal and navigate to the project folder.
2. Navigate to the Simulation folder.
3. Start the run by typing
```bash
     ./runMD.sh <input settings> 
```

4. When the simulation is done in the terminal, all the necessary simulation files will be stored in the _DONE folder.
5. To use the postproccesing module on the raw data generated type
```bash
     ./postprocessing.sh <curent simulation folder> 
```
6. This will generate a csv file containing the calculated quantities.
```csv
T: 303.200984
E_tot: 29.576340
E_kin: 28.570783
E_pot: 1.005557
V: 8609.787683
D: 0.000000
Cv: 295.860007
B: 113.154165
G: 82.268309
E: 198.659941
T_D: 444.161343
```



## Setting.json
The settings file consists of the following inputs:
```json
{
  "Physical_environment" : {"Temperature": 300, "Pressure": 0 },
  "Simulations_config" : {"Timestep": 1,  "Number_of_steps": 10000,  "PBC" : [false, true, false],
  "Potential": { "Kind": "LJ", "Parameters": {"Material": "Ar","epsilon_eV": 0.0103, "sigma": 3.405, "RC": 2.5, "RO": 0.9}},
  "Supercells": [9,9,9], "Tdamp": 100, "Pdamp": 1000, "Friction": 0.05
  },

  "Input_structure" : "../AtomicStructure/test.cif",
  
  "Find_equilibrium" : true,
  "Compute_quantities" : ["E_coh", "B", "Lat_const", "CVT", "Debye", "MSD", "L_crit"]
}

```
### Parameters:
- Physical environment
  - Temperature -  As a positive float (K) * 
  - Pressure -  As a float (Pa) 

  - Timestep - Time step as a float (fs) *
  - Number_of_steps -  As an integer; Total number of steps as an integer *
  - PBC - As a 1x3 list of booleans; Periodic boundary conditions in x,y,z direction
- Potential
  - Kind - As a string;  Which interatomic potential; LJ (Lennard Jones), EMT, MACE *
  - Parameters - Input parameters for Lennard Jones potential incl. Material, epsilon, sigma
- Supercells 
    - Supercell size - As a 1x3 list of integers; If you want a supercell the value 3 gives a 3 x 3 x 3 supercell. For a single unit cell use 0. *  
    - Friction -As a float;  Factor to describe the interaction with the thermostat/heat bath
- Input structure - As a path; The cif/POSCAR file containing the atomic configuration

- Find equilibrium - As a boolean; 
- Compute quantities; As a list of strings; Containing all quantities that should be calculated

For all calculations the settings with (*) are required. For NPT, pressure and compressibility are also required. For NVT, friction is required. If a setting is missing.

### Terminal inputs

positional arguments:
  input_settings                           Path to settings file

options:
  -h, --help                               show this help message and exit
  -PE, --Physical_environment              Dict of temperature and pressure
  -SC, --Simulations_config                Simulation-specific settings
  -IS, --Input_structure                   Initial structure of material to be simulated
  -CQ, --Compute_quantities                List of abbreviations for quantities to compute
  -E, --Ensemble <ENSEMBLE>                Ensemble (NVE, NPT, NVT)
  -T, --Temperature <TEMPERATURE>          Temperature in K
  -P, --Pressure <PRESSURE>                Pressure in Pa
  -PBC <PBC>                               PBC in each direction as a list of booleans
  -POT, --Potential <POTENTIAL>            Potential as a string (EMT, LJ, MACE)
  -TS, --Timestep <TIMESTEP>               Timestep as a float (fs)
  -C, --Compressibility <TIMESTEP>         Compressibility as a float (GPA)
  -µ, --Friction <FRICTION>                Friction coefficent as a float (For NVT)
  -TD, --Tdamp <TDAMP>                     Tdamp as a float (For NPT)
  -PD, --Pdamp <PDAMP>                     Pdamp as a float (For NPT)
  -S, --Supercells <SUPERCELL>             Repetition of input cell e.g [3,3,3], use [1,1,1] for
                                           only unit cell
  -O, --Output_file <PATH>                 Path to where the output file will be written
  -N, --Number_of_steps <NUMBER_OF_STEPS>  Total number of timesteps as an integer
  --debug                                  Debug
  -FE, --Find_equilibrium BOOL             Bool of whether to find equilibrium or not

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

## cif
This file includes the atomic configuration. Example below.

````text
# cif file generated by httk
data_728MSRAF1H
_cell_length_a 10.64217694
_cell_length_b 10.64217694
_cell_length_c 18.61020000
_cell_angle_alpha 90.00000000
_cell_angle_beta 90.00000000
_cell_angle_gamma 120.00000000
_symmetry_space_group_name_hall 'P 1'
_symmetry_space_group_name_h-m 'P 1'
_symmetry_Int_Tables_number 1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_symmetry_multiplicity
_atom_site_Wyckoff_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 1 a 0.00000000 0.16666667 0.50000000 1.00000000 
Mo1 Mo 1 a 0.66666667 0.66666667 0.50000000 1.00000000 
Mo2 Mo 1 a 0.66666667 0.33333333 0.50000000 1.00000000 
Mo3 Mo 1 a 0.66666667 0.00000000 0.50000000 1.00000000 
Mo4 Mo 1 a 0.33333333 0.66666667 0.50000000 1.00000000 
Mo5 Mo 1 a 0.33333333 0.33333333 0.50000000 1.00000000 
Mo6 Mo 1 a 0.33333333 0.00000000 0.50000000 1.00000000 
Mo7 Mo 1 a 0.00000000 0.00000000 0.50000000 1.00000000 
Mo8 Mo 1 a 0.00000000 0.33333333 0.50000000 1.00000000 
Mo9 Mo 1 a 1.00000000 0.66666667 0.50000000 1.00000000 
Te1 Te 1 a 0.22222222 0.77777778 0.40300480 1.00000000 
Te2 Te 1 a 0.88888889 0.77777778 0.59699520 1.00000000 
Te3 Te 1 a 0.22222222 0.11111111 0.59699520 1.00000000 
Te4 Te 1 a 0.88888889 0.44444444 0.40300480 1.00000000 
Te5 Te 1 a 0.88888889 0.44444444 0.59699520 1.00000000 
Te6 Te 1 a 0.22222222 0.11111111 0.40300480 1.00000000 
Te7 Te 1 a 0.88888889 0.11111111 0.40300480 1.00000000 
Te8 Te 1 a 0.88888889 0.11111111 0.59699520 1.00000000 
Te9 Te 1 a 0.55555556 0.77777778 0.40300480 1.00000000 
Te10 Te 1 a 0.55555556 0.77777778 0.59699520 1.00000000 
Te11 Te 1 a 0.22222222 0.44444444 0.59699520 1.00000000 
Te12 Te 1 a 0.55555556 0.44444444 0.40300480 1.00000000 
Te13 Te 1 a 0.88888889 0.77777778 0.40300480 1.00000000 
Te14 Te 1 a 0.22222222 0.44444444 0.40300480 1.00000000 
Te15 Te 1 a 0.55555556 0.11111111 0.40300480 1.00000000 
Te16 Te 1 a 0.55555556 0.11111111 0.59699520 1.00000000 
Te17 Te 1 a 0.22222222 0.77777778 0.59699520 1.00000000 
Te18 Te 1 a 0.55555556 0.44444444 0.59699520 1.00000000 

````

## Further details 
For more details read the technical documentation or
[read our sphinx documentation.](https://isacks-co.github.io/Isacks-Co/)
