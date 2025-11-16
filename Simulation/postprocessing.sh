#!/usr/bin/env bash


folder="$1" 
cd $folder/Outputfiles
python3 ../../../SourceCode/PostProcessing.py "$PWD"



