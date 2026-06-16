import json
import httk.db
from classes import DefectHull, DefectInfo
from abad_classes import AbadParameters
from DBClasses import MDAbadParameters
import matplotlib.pyplot as plt
import numpy as np


def scatter_with_identity_line(x, y, data):
    """
    Creates a scatter plot with the reference line x = y.

    Parameters:
        x (list or array-like): List of x values
        y (list or array-like): List of y values
        title (str): Plot title
    """
    # Convert to numpy arrays
    x = np.array(x)
    y = np.array(y)
    MAE = np.mean(np.abs(x - y))
    RMSE = np.sqrt(np.mean((x - y) ** 2))

    # Check same length
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")

    # Determine limits for identity line
    min_val = min(np.min(x), np.min(y))
    max_val = max(np.max(x), np.max(y))

    # Plot
    plt.figure()
    plt.scatter(
        data["int"]["DFT"],
        data["int"]["MACE"],
        s=5,
        alpha=0.3,
        label="int",
        color="tab:orange",
    )
    plt.scatter(
        data["ads"]["DFT"],
        data["ads"]["MACE"],
        s=5,
        alpha=0.3,
        label="ads",
        color="tab:blue",
    )
    plt.plot(
        [min_val, max_val],
        [min_val, max_val],
        color="grey",
        linestyle="--",
    )
    plt.text(
        0.01,
        0.99,
        f"MAE = {MAE:.3f}",
        transform=plt.gca().transAxes,
        verticalalignment="top",
    )
    plt.text(
        0.01,
        0.95,
        f"RMSE = {RMSE:.3f}",
        transform=plt.gca().transAxes,
        verticalalignment="top",
    )
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.title("Defect hull formation energies")
    plt.xlabel("Formation energy DFT (eV)")
    plt.ylabel("Formation energy MACE (eV)")
    plt.legend()
    plt.grid(True)
    plt.gca().set_aspect("equal", adjustable="box")  # makes x=y visually correct

    # plt.savefig(f"../Formation_energy_plots/{title}")
    plt.show()


with open("formation_energies.json", "r") as f:
    MACE_formation_energies = json.load(f)

with open("converged_structures.json", "r") as f:
    converged_DFT_structures = json.load(f)


def check_defect_type(depth, type):
    if (abs(depth) < 1 and type == "int") or (abs(depth) >= 1 and type == "ads"):
        return True
    else:
        return False


backend = httk.db.backend.Sqlite("../../0K_defects.sqlite")
store = httk.db.store.SqlStore(backend)

search = store.searcher()

defect_hull = search.variable(DefectHull)
DFT_abad = search.variable(AbadParameters)
MACE_abad = search.variable(MDAbadParameters)
info = search.variable(DefectInfo)
search.add(DFT_abad.expansion_factor < 2)
search.add(MACE_abad.expansion_factor < 2)
search.add(defect_hull.defect_key == DFT_abad.key)
search.add(defect_hull.defect_key == MACE_abad.key)
search.add(defect_hull.defect_key == info.key)
search.output(defect_hull, "hull")
search.output(DFT_abad, "DFT_abad")
search.output(MACE_abad, "MACE_abad")
search.output(info, "info")
data = {"int": {"MACE": [], "DFT": []}, "ads": {"MACE": [], "DFT": []}}
MACE_energies = []
DFT_energies = []
i = 0
for match, _ in search:
    formation_energy_DFT = match[0].formation_energies[0]
    formation_energy_MACE = MACE_formation_energies[str(match[0].defect_key)][
        "Formation_energy"
    ]
    depth_DFT = match[1].depth
    depth_MACE = match[2].depth

    defect_type = match[3].defect_type

    if (
        -10 < formation_energy_DFT < 30
        and check_defect_type(depth_DFT, defect_type)
        and defect_type == "int"
    ) or (
        -10 < formation_energy_DFT < 50
        and check_defect_type(depth_MACE, defect_type)
        and defect_type == "ads"
    ):
        print(match[3].host_name)
        data[defect_type]["MACE"].append(formation_energy_MACE)
        data[defect_type]["DFT"].append(formation_energy_DFT)
        DFT_energies.append(formation_energy_DFT)
        MACE_energies.append(formation_energy_MACE)
        diff = abs(formation_energy_DFT - formation_energy_MACE)
        if diff > 1.1:
            i += 1
print(i)
# print(MACE_energies, DFT_energies)
print(len(DFT_energies))
scatter_with_identity_line(DFT_energies, MACE_energies, data)
