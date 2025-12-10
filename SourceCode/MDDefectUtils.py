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

from abad_classes import Delta
from classes import DefectCell, ScreenResult

from DBClasses import MDDelta, MDScreenResult


def getAllDefectKeys(store):
    """
    Function to get all defect keys from DefectCell.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
    Returns:
        all_Keys (list) : A list of all the defect keys from the DefectCell table.
    """
    search = store.searcher()
    search_defectcell = search.variable(DefectCell)

    search.output(search_defectcell.host_name, "host")
    search.output(search_defectcell.defect_name, "defect")
    search.output(search_defectcell.key, "defect_key")

    all_keys = set()
    for match, header in search:
        host_name = match[0]
        defect_name = match[1]
        key = match[2]

        all_keys.add((host_name, defect_name, key))

    return all_keys


def getMDKeys(store):
    """
    Function to get all keys from MDScreenResult table.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
    Returns:
        all_md_keys (list) : A list of all the defect keys from the MDScreenResult table.
    """
    search = store.searcher()
    search_MD_cell = search.variable(MDScreenResult)

    search.output(search_MD_cell.defect_key, "defect_key")

    all_md_keys = set()
    for match, header in search:
        key = match[0]

        all_md_keys.add(key)

    return all_md_keys


def getDelta(store):
    """
    Function to get all delta values from the Delta table.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
    Returns:
        delta (list) : A list containing all the host and defect names from the Delta table.
    """
    search = store.searcher()
    search_delta = search.variable(Delta)
    search.output(search_delta.host, "host")
    search.output(search_delta.dopant, "dopant")

    all_delta = set()

    for match, header in search:
        host = match[0]
        defect = match[1]
        all_delta.add((host, defect))

    return all_delta


def getMDDelta(store):
    """
    Function to get all delta values from the MDDelta table.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
    Returns:
        MDDelta (list) : A list of tuples containing all the finished delta calculations delta_host and delta_defect.
    """
    search = store.searcher()
    search_MD_delta = search.variable(MDDelta)
    search.output(search_MD_delta.host, "host")
    search.output(search_MD_delta.dopant, "dopant")

    all_MD_delta = set()

    for match, header in search:
        delta_host = match[0]
        delta_dopant = match[1]
        all_MD_delta.add((delta_host, delta_dopant))

    return all_MD_delta


def notDoneMDRuns(store):
    """
    Function to compare all keys from DefectCell and MDScreenResult to find out which simulations that have not been
    completed.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
    Returns:
        missing_keys (set) : A set of keys belonging to simulations not completed.
    """
    defect_keys = getAllDefectKeys(store)
    md_keys = getMDKeys(store)

    missing_keys = set(info for info in defect_keys if not info[2] in md_keys)

    return missing_keys


def notDoneDelta(store):
    """
    Function to compare all delta values between DefectCell and MDScreenResult to find out which simulations that
    have not been completed.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
    Returns:
        missing_keys (set) : A set of keys that have not been completed.
    """
    delta = getDelta(store)
    MD_delta = getMDDelta(store)

    missing_delta = delta - MD_delta

    return missing_delta
