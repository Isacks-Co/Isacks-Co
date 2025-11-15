#!/usr/bin/env bash

file1="$1"
file2="$2"

# Extract extension (everything after the last dot)
ext="${file1##*.}"

case "$ext" in
    
    
    json)
        
        python3 ../SourceCode/scriptUtils/setup.py "$file2" "$file1"
        ;;
    *)
        python3 ../SourceCode/scriptUtils/setup.py "$file1" "$file2"
        ;;
esac
