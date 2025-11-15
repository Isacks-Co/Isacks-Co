import sys
import os 
import re
if __name__ == "__main__":
   
  
    

    outputfile = sys.argv[1]

    # Get the current working directory (the folder you want to rename)
    dir_path = os.getcwd()

    # Read the first line from the file
    with open(outputfile, "r") as f:
        first_line = f.readline().strip()

    # Sanitize the string to remove invalid characters for folder names
    safe_name = re.sub(r'[\/:*?"<>|]', "_", first_line)

    # Build the full new path
    new_path = os.path.join(os.path.dirname(dir_path), f"{safe_name}_DONE")

    # Rename the folder
    os.rename(dir_path, new_path)
