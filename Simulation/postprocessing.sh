#!/usr/bin/env bash


folder="$1" 
cd $folder/OutputFiles
python3 ../../../SourceCode/PostProcessing.py "$PWD" "../SetupFiles/settings.json"



