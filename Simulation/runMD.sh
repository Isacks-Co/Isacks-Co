#!/bin/bash

file1="$1"
file2="$2"
# Get current directory
current_dir=$(pwd)

bash setup.sh $file1 $file2

cd currentSimulation/OutputFiles

python3 ../../../SourceCode/MolecularDynamics.py ../SetupFiles/POSCAR ../SetupFiles/settings.json

python3 ../../../SourceCode/PostProcessing.py "$PWD"

cd ..

python3 ../../SourceCode/scriptUtils/renameFolder.py OutputFiles/Quantities.csv