#!/usr/bin/env bash


# Get current directory
current_dir=$(pwd)

# Loop over all directories in the current directory
for dir in "$current_dir"/*/; do
    if [ -d "$dir" ]; then
        echo "Deleting directory: $dir"
        rm -rf "$dir"
    fi
done

echo "All directories in $current_dir have been deleted."
