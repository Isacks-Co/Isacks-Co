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
from classes import *

from MDDefectUtils import notDoneMDRuns, notDoneDelta

INF = 10000000

backend = httk.db.backend.Sqlite('../../defect/defects.sqlite')
store = httk.db.store.SqlStore(backend)
search = store.searcher()

search_hosts = search.variable(DefectInfo)
search.output(search_hosts, "hosts")

host_list = []
dopant_list = []

# Find all hosts materials and dopant in database
for match in search:
    host = str(match[0][0].host_name)
    dopant = re.match(r"^([^_]+)", match[0][0].defect_name).group(1)

    if host not in host_list:
        host_list.append(host)

    if dopant not in dopant_list:
        dopant_list.append(dopant)

missing_runs = notDoneMDRuns(store)
missing_delta = notDoneDelta(store)

# Loop over all hosts in the database
for host in host_list:
    search = store.searcher()
    search_defectinfo = search.variable(DefectInfo)
    search_screenresult = search.variable(ScreenResult)
    search_abad = search.variable(AbadParameters)

    search.add(search_defectinfo.key == search_screenresult.defect_key)
    search.add(search_defectinfo.key == search_abad.key)
    search.add(search_defectinfo.host_name == host)

    search.output(search_defectinfo, 'defectinfo')
    search.output(search_screenresult, 'screenresult')
    search.output(search_abad, 'abad')

    # Loop over a specific host and dopant combination
    for dopant in dopant_list:

        # Check if we are missing any simulation for combinaton of host and dopant
        if any((missing[0] == host and re.match(r"^([^_]+)", missing[1]).group(1) == dopant) for missing in
               missing_runs):
            print(f"{host} : {dopant} are missing simulations. No Delta calculated")
            continue

        # Check if we already have a stored value in the database
        if (host[:-4], dopant) not in missing_delta:
            continue

        adatom_lowest, interstital_lowest = INF, INF
        adatom_dopant, interstital_dopant = "", ""
        for match in search:
            defect_info = match[0][0]
            screenresult = match[0][1]
            abad = match[0][2]

            if re.match(r"^([^_]+)", defect_info.defect_name).group(1) != dopant:
                continue

            # Skip values from simulation that exploded
            if abad.expansion_factor > 1.8:
                continue

            if re.match(r"^.*_ads\d+$", defect_info.defect_name) and screenresult.total_energy_coarse < adatom_lowest:
                adatom_lowest = screenresult.total_energy_coarse
                adatom_dopant = defect_info.defect_name

            if re.match(r"^.*_int\d+$",
                        defect_info.defect_name) and screenresult.total_energy_coarse < interstital_lowest:
                interstital_lowest = screenresult.total_energy_coarse
                interstital_dopant = defect_info.defect_name

        # Make sure that we found some values
        if adatom_lowest == INF or interstital_lowest == INF:
            print(f"{host} : {dopant} No valid results found in the database, delta not calculated")
            continue

        # print(f"For the host : {defect_info.host_name}, the lowest adatom, respectively interstital, energy is "
        #       f"\n{adatom_lowest} eV for {adatom_dopant} and {interstital_lowest} eV for {interstital_dopant} eV for {interstital_dopant}")

        # TODO Fix upload logisics
        delta = interstital_lowest - adatom_lowest

        print(
            f"Delta for {host} - {dopant} : {delta} eV, with material {interstital_dopant if delta < 0 else adatom_dopant}.")
