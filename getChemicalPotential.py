# from mace.calculators import MACECalculator

from ASEWrappers import LangevinIntegrator, MACEPotential, AtomicStructure
from MDClasses import EquilibriumRun
import SimulationInput
from pathlib import Path

"""
def getChemicalPotential(atoms):
    # bulk_material = bulk(element)
    atoms.calc = MACECalculator(
        model_paths="mace-mpa-0-medium.model",
        device="cpu",
        default_dtype="float64",
        head="default",
    )
    # bulk_material = AtomicStructure(bulk_material)
    mu = atoms.get_potential_energy() / len(atoms)
    return mu"""


def relax_defects(paths):
    num_steps = 500
    temp_k = 0
    time_steps = 1.0
    pot = MACEPotential("mace-mpa-0-medium.model")

    data = {}
    for i, p in enumerate(paths):
        print(i)
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
        # print(p.split("/")[1][:-7])

        data[p.split("/")[1][:-7]] = E_post / len(atomic_structure)
        if abs(E_post - E_pre) > 1e-1:
            print(E_pre, E_post, p.split("/")[1][:-7])

    return data


from mp_api.client import MPRester

API_KEY = "0XFSEY8CIPcmgtadsHJfs12PkKDuCkOI"


def get_most_stable_unary_ase(element_symbol):
    """
    Returns the most stable unary material for a given element as an ASE Atoms object.
    Returns None if no stable unary material exists.
    """
    with MPRester(API_KEY) as mpr:
        # Query stable unary materials
        results = mpr.materials.summary.search(
            elements=[element_symbol],
            num_elements=1,
            # energy_above_hull=(0, None),
            fields=["material_id", "formula_pretty", "energy_above_hull"],
        )

        if len(results) == 0:
            return None
        # Sort to get the most stable polymorph
        most_stable = sorted(results, key=lambda x: x.energy_above_hull)[0]
        # print(most_stable.formula_pretty)
        # Get pymatgen Structure
        structure = mpr.materials.get_structure_by_material_id(most_stable.material_id)
        structure.to(filename=f"defects/{element_symbol}.poscar", fmt="poscar")
        # Convert to ASE Atoms
        # atoms = AseAtomsAdaptor.get_atoms(structure)

        # return atoms


def loop_files(folder_path):
    paths = []
    for file_path in Path(folder_path).iterdir():
        if file_path.is_file():
            print(file_path)
            paths.append(str(file_path))

    return paths


if __name__ == "__main__":
    # paths = loop_files("missing_defects")
    r = relax_defects(["missing_defects/W.poscar"])
    print(r)
    quit()

    """paths = loop_files("defects")
    with open("relaxed_chem_pot.json", "r") as f:
        data = json.load(f)
    pot = MACEPotential("mace-mpa-0-medium.model")
    res = {}
    for p in paths:
        atomic_structure = AtomicStructure.fromFile(path=p, potential=pot)
        N = len(atomic_structure)
        print(N)
        res[p.split("/")[1][:-7]] = data[p.split("/")[1][:-7]] / N
    with open("relaxed_chem_pot2.json", "w") as f:
        json.dump(res, f)"""
