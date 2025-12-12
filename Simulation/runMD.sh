#!/usr/bin/env bash


file1="$1"
shift
args="$@"
# Get current directory
current_dir=$(pwd)

bash setup.sh $file1

cd currentSimulation/OutputFiles
touch sampledata.txt

python3 ../../../SourceCode/MolecularDynamics.py ../SetupFiles/settings.json $args



cd ../

python3 ../../SourceCode/scriptUtils/renameFolder.py OutputFiles/sampledata.txt
