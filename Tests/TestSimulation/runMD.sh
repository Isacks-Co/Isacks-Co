#!/usr/bin/env bash


file1="$1"
shift
args="$@"
# Get current directory
current_dir=$(pwd)

bash setup.sh $file1

cd currentSimulation/OutputFiles

python3 ../../../../SourceCode/MolecularDynamics.py ../SetupFiles/settings.json $args

if [ $? -eq 0 ]; then
    
    cd ../..
    bash postprocessing.sh currentSimulation
    cd currentSimulation/OutputFiles
    cd ..
    touch OutputFiles/sampledata.txt
    python3 ../../../SourceCode/scriptUtils/renameFolder.py OutputFiles/sampledata.txt
    
else
    exit 1 # Exit with a generic failure code
fi

