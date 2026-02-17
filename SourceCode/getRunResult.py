import json
import bz2
import pymongo
from MDStoreUtils import saveMDScreenResult, saveMDAbadParameters, CommitAndClose, saveMDQuantities
import sys
import os
import httk.db
from ase.io import read
from optimade.adapters.structures.ase import from_ase_atoms
from optimade.models.structures import StructureResource
from datetime import datetime, timezone

client = pymongo.MongoClient("mongodb://localhost:27017/")
db = client["molecular_dynamics_db"]
collection = db["structures"]

root = sys.argv[1]          #Path to top folder, test: ../../test
db_path = sys.argv[2]       # 🥀💀 for me, personally '../../defect/defects.sqlite'

backend = httk.db.backend.Sqlite(db_path)
store = httk.db.store.SqlStore(backend)

current_time = datetime.now(timezone.utc).isoformat()

optimade_list = []
leaf_dirs = []
for current_dir, subdirs, files in os.walk(root):
    if not subdirs:
        leaf_dirs.append(current_dir)

for d in leaf_dirs:
    key, energy = 0,0
    atoms = None
    optimade_attributes = None
    for file in os.listdir(d):
        fullpath = os.path.join(d, file)
        if file.endswith(".json.bz2"):

            with bz2.open(fullpath) as json_file:
                data = json.load(json_file)

                if "MDScreenResult" in data:
                    info = data["MDScreenResult"]
                    key = info["key"]
                    energy = info["energy"]
                    saveMDScreenResult(
                        store =store,
                        defect_key=key,
                        total_energy=energy)

                if "MDAbadParameters" in data:
                    info = data["MDAbadParameters"]
                    saveMDAbadParameters(
                        store=store,
                        defect_key=info["key"],
                        depth=info["depth"],
                        expansion_factor=info["expansion_factor"],
                        defect_index=info["defect_index"],
                        lattice_constant=info["lattice_constant"],
                        convergence_criterium=info["convergence_criterium"],
                        number_of_steps=info["number_of_steps"],
                        time=info["time"],
                    )

                if "MDQuantities" in data:
                    info = data["MDQuantities"]
                    saveMDQuantities(
                            store = store,
                            key=info["key"],
                            temperature=info["temperature"],
                            energy_pot=info["energy_pot"],
                            energy_kinetic=info["energy_kinetic"],
                            msd=info["msd"],
                            internal_pressure=info["internal_pressure"],
                        )

        if file.startswith("post_structure"):
            atoms = read(fullpath)
            # Remove extra information that cif file provides. Gets stored as self made attributes otherwise.
            for key in list(atoms.info):
                del atoms.info[key]
            # Stores the data from atoms object as a ResourceStructureAttributes optimade object

    if (atoms is None) or (key == 0 or energy == 0):
        print(f"The directory {d} are missing results for creating OPTIMADE entry into DB")
        continue

    optimade_attributes = from_ase_atoms(atoms=atoms)
    optimade_attributes = optimade_attributes.model_dump(mode="json")

    # OPTIMADE requires structure_features to be a list, never null
    if optimade_attributes.get("structure_features") is None:
        optimade_attributes["structure_features"] = []                      # Kept getting errors due to being written as Null
    optimade_attributes["last_modified"] = current_time
    optimade_attributes["_MD_total_energy"] = energy

    structure_resource = StructureResource(id=str(key), type="structures", attributes=optimade_attributes)
    optimade_list.append(structure_resource.model_dump())       #optimade_list.append(structure_resource.model_dump(exclude_unset=True))

CommitAndClose(backend)
result = collection.insert_many(optimade_list, ordered=False)
