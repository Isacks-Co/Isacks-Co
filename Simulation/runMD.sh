#!/usr/bin/env bash


file1="$1"
flag="$2"
# Get current directory
current_dir=$(pwd)

bash setup.sh $file1

cd currentSimulation/OutputFiles
touch sampledata.txt

python3 ../../../SourceCode/MolecularDynamics.py ../SetupFiles/settings.json






if [[ $flag == "-P" ]]; then
    cd ../..
    bash postprocessing.sh currentSimulation
    cd currentSimulation/OutputFiles
fi

cd ../

python3 ../../SourceCode/scriptUtils/renameFolder.py OutputFiles/sampledata.txt