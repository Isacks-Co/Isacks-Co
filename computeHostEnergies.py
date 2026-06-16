"""db_path = "../../0K_defects.sqlite"
backend = httk.db.backend.Sqlite(db_path)
store = httk.db.store.SqlStore(backend)

search = store.searcher()

data = {}
hosts = search.variable(HostSuperCell)
search.output(hosts, "host")

for match, _ in search:
    host = match[0].host_supercell
    # print(structure_to_ase_atoms(host))
    host.io.save(f"hosts/{match[0].name()}.cif")
"""

from pathlib import Path

from ASEWrappers import LangevinIntegrator, MACEPotential, AtomicStructure
from MDClasses import EquilibriumRun
import SimulationInput


def loop_files(folder_path):
    paths = []
    for file_path in Path(folder_path).iterdir():
        if file_path.is_file():
            print(file_path)
            paths.append(str(file_path))

    return paths


def relax_hosts(paths):
    num_steps = 500
    temp_k = 0
    time_steps = 1.0
    pot = MACEPotential("mace-mpa-0-medium.model")

    data = {}
    for p in paths:
        atomic_structure = AtomicStructure.fromFile(path=p, potential=pot)
        E_pre = atomic_structure.potential_energy
        print(E_pre)
        lang_int = LangevinIntegrator(time_steps, temp_k, 0.05)
        settings = SimulationInput.SimulationSettings(
            num_steps, pot, lang_int, "NVT", False
        )
        equil_MD = EquilibriumRun(settings=settings)
        equil_structure = equil_MD.run(
            atomic_structure,
            settings.num_steps,
            store_traj=False,
            check_conv=True,
            check_expansion=False,
        )

        E_post = equil_structure.potential_energy
        print(E_post)
        data[p.split("/")[1][:-8]] = E_post
        if abs(E_post - E_pre) > 1e-1:
            print(E_pre, E_post, p.split("/")[1][:-8])
        # print(data)

    return data


import json

if __name__ == "__main__":
    paths = loop_files("host_super")

    relaxed_energies = relax_hosts(paths)

    with open("relaxed_host_energies.json", "w") as f:
        json.dump(relaxed_energies, f, indent=4)
