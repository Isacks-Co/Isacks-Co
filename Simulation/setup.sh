#!/usr/bin/env bash

settings_file="$1"

# Extract extension (everything after the last dot)
ext="${settings_file##*.}"

case "$ext" in
    
    
    json)
        
        python3 ../SourceCode/scriptUtils/setup.py "$settings_file"
        ;;
    *)
        python3 ../SourceCode/scriptUtils/setup.py "$settings_file"
        ;;
esac
