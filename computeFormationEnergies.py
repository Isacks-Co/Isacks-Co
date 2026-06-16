import httk.db
import json
from DBClasses import MDScreenResult
from classes import DefectInfo

with open("relaxed_host_energies.json", "r") as f:
    host_energies = json.load(f)
with open("relaxed_chem_pot2.json", "r") as f:
    chem_pot = json.load(f)

db_path = "../../0K_defects.sqlite"
backend = httk.db.backend.Sqlite(db_path)
store = httk.db.store.SqlStore(backend)

search = store.searcher()

screen = search.variable(MDScreenResult)
info = search.variable(DefectInfo)

search.add(screen.defect_key == info.key)
search.output(screen, "s")
search.output(info, "c")
res = {}
missing = []
for match, _ in search:
    defect_bulk_energy = match[0].total_energy_coarse
    # print(match[0].key)
    host = match[1].host_name[:-4]
    defect = match[1].defect_name
    # print(host, defect)

    formation_energy = defect_bulk_energy - host_energies[host] - chem_pot[defect[:-5]]
    res[str(match[0].defect_key)] = {
        "host": host,
        "defect": defect,
        "defect_type": match[1].defect_type,
        "Formation_energy": formation_energy,
    }

print(missing)

with open("formation_energies.json", "w") as f:
    json.dump(res, f)
