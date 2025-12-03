import os
import sys
import json
from MDStoreUtils import saveMDScreenResult, saveMDAbadParameters, CommitAndClose
import httk.db

root = sys.argv[1] #Path to top folder, test: ../../test
db_path = sys.argv[2] # 🥀💀 for me, personally '../../defect/defects.sqlite'

backend = httk.db.backend.Sqlite(db_path)
store = httk.db.store.SqlStore(backend)


leaf_dirs = []
for current_dir, subdirs, files in os.walk(root):
    if not subdirs:
        leaf_dirs.append(current_dir)

for d in leaf_dirs:
    for file in os.listdir(d):
        fullpath = os.path.join(d, file)
        if file.endswith(".json"):
            #print(fullpath)

            with open(fullpath) as json_file:
                data = json.load(json_file)
                #print(data)

                if "MDScreenResult" in data:
                    info = data["MDScreenResult"]
                    saveMDScreenResult(
                        store =store,
                        defect_key=info["key"],
                        total_energy=info["total_energy_coarse"])

                if "MDAbadParameters" in data:
                    info = data["MDAbadParameters"]
                    saveMDAbadParameters(
                        store=store,
                        defect_key=info["key"],
                        depth=info["depth"],
                        expansion_factor=info["expansion_factor"],
                        defect_index=info["defect_index"])

