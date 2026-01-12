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

import re

import httk.db
from abad_classes import *
from classes import DefectInfo, ScreenResult

from DBClasses import MDDelta, MDScreenResult, MDAbadParameters
from MDDefectUtils import notDoneMDRuns, notDoneDelta
from MDStoreUtils import CommitAndClose


def getDopant(defect_name):
    """
    Extract dopant identifier from a defect name string.

    The dopant is assumed to be the prefix of the defect name before the first
    underscore.

    Parameters
    ----------
    defect_name : str
        Defect name string, e.g. ``"C_ads0"`` or ``"C_int2"``.

    Returns
    -------
    str
        Dopant identifier (prefix before the first underscore).
    """
    return re.match(r"^([^_]+)", defect_name).group(1)


INF = 10000000
# Open SQLite database and create store/searcher interface
backend = httk.db.backend.Sqlite('../../2D_defects.sqlite')
store = httk.db.store.SqlStore(backend)
search = store.searcher()

search_hosts = search.variable(DefectInfo)
search.output(search_hosts, "hosts")
# Collect all unique host materials and dopant identifiers present in DefectInfo
host_list = []
dopant_list = []

# Find all hosts materials and dopant in database
for match in search:
    info = match[0][0]
    host = str(info.host_name)
    dopant = getDopant(info.defect_name)

    if host not in host_list:
        host_list.append(host)

    if dopant not in dopant_list:
        dopant_list.append(dopant)
# Determine which MD runs and deltas are missing so we can skip incomplete pairs
missing_runs = notDoneMDRuns(store)
missing_delta = notDoneDelta(store)

# Loop over all hosts in the database
for host in host_list:
    search = store.searcher()
    search_defectinfo = search.variable(DefectInfo)
    search_screenresult = search.variable(MDScreenResult)
    search_abad = search.variable(MDAbadParameters)

    search.add(search_defectinfo.key == search_screenresult.defect_key)
    search.add(search_defectinfo.key == search_abad.key)
    search.add(search_defectinfo.host_name == host)

    search.output(search_defectinfo, 'defectinfo')
    search.output(search_screenresult, 'screenresult')
    search.output(search_abad, 'abad')

    # collect everything once in list, instead of doing match every single time which is 1 query per host PER defect
    combined_results = []
    for match in search:
        defect_info = match[0][0]
        screenresult = match[0][1]
        abad = match[0][2]
        combined_results.append((defect_info, screenresult, abad))

    # Loop over a specific host and dopant combination
    for dopant in dopant_list:

        # Check if we are missing any simulation for combinaton of host and dopant
        if any(missing_host == host and getDopant(missing_defect) == dopant
               for missing_host, missing_defect, key in missing_runs):
            print(f"{host} : {dopant} are missing simulations. No Delta calculated")
            continue

        # Check if we already have a stored value in the database
        if (host[:-4], dopant) not in missing_delta:
            continue

        adatom_lowest, interstital_lowest = INF, INF
        adatom_defect, interstital_defect = "", ""
        adatom_key, interstital_key = "", ""

        for defect_info, screenresult, abad in combined_results:

            if getDopant(defect_info.defect_name) != dopant:
                continue

            # Skip values from simulation that exploded
            if abad.expansion_factor > 1.8:
                continue

            if re.match(r"^.*_ads\d+$",
                        defect_info.defect_name) and screenresult.total_energy_coarse < adatom_lowest:
                adatom_lowest = screenresult.total_energy_coarse
                adatom_defect = defect_info.defect_name
                adatom_key = defect_info.key

            if re.match(r"^.*_int\d+$",
                        defect_info.defect_name) and screenresult.total_energy_coarse < interstital_lowest:
                interstital_lowest = screenresult.total_energy_coarse
                interstital_defect = defect_info.defect_name
                interstital_key = defect_info.key

        # Make sure that we found some values
        if adatom_lowest == INF or interstital_lowest == INF:
            print(f"{host} : {dopant} No valid results found in the database, delta not calculated")
            continue

        delta = interstital_lowest - adatom_lowest

        # print(
        #     f"Delta for {host} - {dopant} : {delta} eV, with material {interstital_defect if delta < 0 else adatom_defect}.")

        store.save(MDDelta(
            host=host,
            dopant=dopant,
            defect=interstital_defect if delta < 0 else adatom_defect,
            key=interstital_key if delta < 0 else adatom_key,
            delta=delta,
        ))

CommitAndClose(backend)
