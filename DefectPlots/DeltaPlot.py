# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


#!/usr/bin/env python3
from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from ase.db import connect

from httk.db.backend import Sqlite
from httk.db.store import SqlStore
import sys
sys.path.append("../SourceCode")

from DBClasses import MDDelta, MDAbadParameters



# ============================================================
# CONFIG: paths to current databases
# ============================================================
DEFECTS_DB_PATH = "../../2D_defects.sqlite"

# ============================================================
# 1. Periodic table groups and period panels (Fig. 3 layout)
# ============================================================
GROUP_TABLE = {
    # Period 1
    "H": 1,  "He": 18,

    # Period 2
    "Li": 1, "Be": 2, "B": 13, "C": 14, "N": 15, "O": 16, "F": 17, "Ne": 18,

    # Period 3
    "Na": 1, "Mg": 2, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,

    # Period 4
    "K": 1,  "Ca": 2,
    "Sc": 3, "Ti": 4, "V": 5,  "Cr": 6,  "Mn": 7,  "Fe": 8,  "Co": 9,  "Ni": 10,
    "Cu": 11, "Zn": 12, "Ga": 13, "Ge": 14, "As": 15, "Se": 16, "Br": 17, "Kr": 18,

    # Period 5
    "Rb": 1, "Sr": 2,
    "Y": 3,  "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10,
    "Ag": 11, "Cd": 12, "In": 13, "Sn": 14, "Sb": 15, "Te": 16, "I": 17, "Xe": 18,

    # Period 6 (s + p + 5d, Ln separate)
    "Cs": 1, "Ba": 2,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10,
    "Au": 11, "Hg": 12, "Tl": 13, "Pb": 14, "Bi": 15, "Po": 16, "At": 17, "Rn": 18,

    # Lanthanides (own "period")
    "La": 4, "Ce": 5, "Pr": 6, "Nd": 7, "Pm": 8, "Sm": 9,
    "Eu": 10, "Gd": 11, "Tb": 12, "Dy": 13, "Ho": 14, "Er": 15,
    "Tm": 16, "Yb": 17, "Lu": 18,
}

PERIOD_PANELS = [
    ["H", "He"],  # period 1
    ["Li", "Be", "B", "C", "N", "O", "F", "Ne"],  # period 2
    ["Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"],  # period 3
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
     "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],  # period 4
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
     "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],  # period 5
    ["Cs", "Ba", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
     "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],  # period 6
    ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
     "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"],  # lanthanides
]


# ============================================================
# 2. CURRENT DATA SOURCE: defects.sqlite + imp2d.db
# ============================================================

def load_deltas_from_defects(defects_db=DEFECTS_DB_PATH):
    """
    CURRENT IMPLEMENTATION:
        Δ(H,X) from defects.sqlite via Delta.delta

    LATER (when you have your own data):
        Replace THIS function with one that returns dopants, hosts, deltas
        from your own source instead of defects.sqlite.
    """
    backend = Sqlite(defects_db)
    store = SqlStore(backend)

    search = store.searcher()
    search_delta = search.variable(MDDelta)
    search_abad = search.variable(MDAbadParameters)
    search.add(search_delta.key == search_abad.key)
    search.output(search_delta, "MDDelta")
    search.output(search_abad, "MDAbadParameters")

    dopants = []
    hosts = []
    delta_vals = []
    lattice_const = []

    for match, header in search:
        delta_obj = match[0]
        abad_obj = match[1]

        dopants.append(delta_obj.dopant)
        hosts.append(delta_obj.host)
        delta_vals.append(float(delta_obj.delta))
        lattice_const.append(float(abad_obj.lattice_constant))

    return np.array(dopants), np.array(hosts), np.array(delta_vals, dtype=float), np.array(lattice_const, dtype=float)


def prepare_delta_plot_data(
    defects_db=DEFECTS_DB_PATH
):
    """
    Combine Δ(H,X) and a_mean(host) into arrays for plotting.

    CURRENT:
      - Δ from defects.sqlite
    """
    dopants_raw, hosts_raw, deltas_raw, lattice_raw = load_deltas_from_defects(defects_db)

    dopants = []
    deltas = []
    a_means = []

    for d, h, val, a in zip(dopants_raw, hosts_raw, deltas_raw, lattice_raw):
        # element must be in GROUP_TABLE (we only plot those)
        if d not in GROUP_TABLE:
            continue

        if a is None:
            continue
        a /= 3

        dopants.append(d)
        deltas.append(val)
        a_means.append(a)

    return (
        np.array(dopants),
        np.array(deltas, dtype=float),
        np.array(a_means, dtype=float),
    )


# ============================================================
# 3. Plot: period panels, with trend line
# ============================================================

def plot_by_period_panels_from_delta(
    filename="DeltaPlot.png",
    defects_db=DEFECTS_DB_PATH
):
    dopants, deltas, a_means = prepare_delta_plot_data(defects_db)

    if deltas.size == 0:
        print("No Δ data found – check filters or DB paths.")
        return

        # Remove outliers: keep only -8 <= Δ <= 14
    delta_mask = (deltas >= -8.0) & (deltas <= 14.0)

    dopants_np = dopants[delta_mask]
    deltas_np  = deltas[delta_mask]
    a_means_np = a_means[delta_mask]

    # group position + jitter
    groups_all = np.array([GROUP_TABLE[d] for d in dopants_np])
    rng = np.random.default_rng(0)
    jitter_all = rng.uniform(-0.15, 0.15, size=len(groups_all))
    x_all = groups_all + jitter_all

    # a_mean color normalization: 2–8 Å, clipped
    norm = Normalize(vmin=2.0, vmax=8.0)
    # Here you can adjust the sizes of the graphs
    n_panels = len(PERIOD_PANELS)
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(10, 2 * n_panels),
        sharex=False,
        sharey=True,

    )

    if n_panels == 1:
        axes = [axes]

    top_axis = True
    for ax, panel in zip(axes, PERIOD_PANELS):
        panel_dopants_for_ticks = [d for d in panel if d in GROUP_TABLE]
        panel_dopants_with_data = [
            d for d in panel_dopants_for_ticks if d in dopants_np
        ]

        # Mask for this panel
        if panel_dopants_with_data:
            mask = np.isin(dopants_np, panel_dopants_with_data)
        else:
            mask = np.zeros_like(dopants_np, dtype=bool)

        sc = None

        if mask.any():
            x_vals = x_all[mask]
            y_vals = deltas_np[mask]
            c_vals_raw = a_means_np[mask]
            c_vals = np.clip(c_vals_raw, 2.0, 8.0)

            sc = ax.scatter(
                x_vals,
                y_vals,
                c=c_vals,
                cmap="viridis",
                norm=norm,
                s=20,
                edgecolors="none",
            )

            # ---- AVERAGE BOXES + TREND LINE ----
            line_x = []
            line_y = []

            for d in panel_dopants_with_data:
                g = GROUP_TABLE[d]
                d_mask = (dopants_np == d)
                vals = deltas_np[d_mask]
                if len(vals) == 0:
                    continue
                mean_delta = float(vals.mean())

                # orange square at the mean
                ax.errorbar(
                    g,
                    mean_delta,
                    fmt="s",
                    mfc="none",
                    mec="orange",
                    ecolor="orange",
                    ms=6,
                    lw=1.0,
                    zorder=3,
                )

                line_x.append(g)
                line_y.append(mean_delta)

            # connect the squares with a trend line
            if line_x:
                # sort by group
                line_pairs = sorted(zip(line_x, line_y), key=lambda p: p[0])

                # build segments where consecutive groups differ by 1
                seg_x = [line_pairs[0][0]]
                seg_y = [line_pairs[0][1]]

                for (g_prev, y_prev), (g_curr, y_curr) in zip(line_pairs, line_pairs[1:]):
                    if g_curr - g_prev == 1:
                        # continue same segment
                        seg_x.append(g_curr)
                        seg_y.append(y_curr)
                    else:
                        # draw old segment (if it has at least 2 points)
                        if len(seg_x) > 1:
                            ax.plot(seg_x, seg_y, color="orange", linewidth=1.0, zorder=2)
                        # start new segment
                        seg_x = [g_curr]
                        seg_y = [y_curr]

                # draw last segment
                if len(seg_x) > 1:
                    ax.plot(seg_x, seg_y, color="orange", linewidth=1.0, zorder=2)


        # Δ = 0 reference
        ax.axhline(0.0, linestyle="--", color="0.4", linewidth=1)

        ax.set_ylim(-8, 14)
        ax.set_yticks([0, 10])
        ax.set_xlim(0.5, 18.5)
        ax.set_ylabel(r"$\Delta$ [eV]")

        # bottom: element symbols
        tick_positions = [GROUP_TABLE[d] for d in panel_dopants_for_ticks]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(panel_dopants_for_ticks)

        # top: group numbers
        if top_axis:
            top = ax.secondary_xaxis("top")
            top.set_xticks(range(1, 19))
            top_axis = False

        # Colorbar if there is data
        if sc is not None:
            cbar = fig.colorbar(
            sc,
            ax=ax,
            location="right",
            fraction=0.02,  # thinner bar
            pad=0.01,       # closer to the axes
        )
        cbar.set_label(r"$a_\mathrm{mean}$ [Å]")
        cbar.set_ticks([4, 6])
        cbar.set_ticklabels(["4", "6"])

    plt.tight_layout(rect=[0.06, 0.04, 0.94, 0.96])
    plt.savefig(filename, dpi=300)
    # plt.show()


if __name__ == "__main__":
    # change the filename, defects database is the database path, which right now is defaulted
    # imp2d_db contains all the data (used to compute)
    plot_by_period_panels_from_delta(
        filename="Deltaplot.png",
        defects_db=DEFECTS_DB_PATH,
    )
