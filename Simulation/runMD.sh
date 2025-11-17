#!/usr/bin/env bash


file1="$1"
file2="$2"
flag="$3"
# Get current directory
current_dir=$(pwd)

bash setup.sh $file1 $file2

cd currentSimulation/OutputFiles

python3 ../../../SourceCode/MolecularDynamics.py ../SetupFiles/POSCAR ../SetupFiles/settings.json

if [[ $flag == "-P" ]]; then
    cd ../..
    bash postprocessing.sh currentSimulation
    cd currentSimulation/Outputfiles
fi

cd ../

python3 ../../SourceCode/scriptUtils/renameFolder.py OutputFiles/sampledata.txt