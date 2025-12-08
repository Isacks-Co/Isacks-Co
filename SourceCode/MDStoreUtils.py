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

from DBClasses import MDDefectCell, MDScreenResult, MDAbadParameters


def saveMDDefectCell(store, host_name, defect_name, defect_key, structure):
    """
    Wrapper to save MDDefectCell object to backend database.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
        host_name (str): The name of the host material.
        defect_name (str): The name of the defect material.
        defect_key (str): The key of the defect material.
        structure (MDDefectCell): the structure of the defect material.
    """
    cell = MDDefectCell(
        host_name=host_name,
        defect_structure=structure,
        defect_name=defect_name,
        key=int(defect_key)
    )
    store.save(cell)


def saveMDScreenResult(store, defect_key, total_energy):
    """
    Wrapper to save MDScreenResult object to backend database.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
        defect_key (str): The key of the defect material.
        total_energy (float): The total energy of the defect material.
    """
    result = MDScreenResult(
        defect_key=int(defect_key),
        total_energy_coarse=total_energy,
    )
    store.save(result)


def saveMDAbadParameters(store, defect_key, depth, expansion_factor, defect_index, lattice_constant):
    """
    Wrapper to save MDAbadParameters object to backend database.
    Args:
        store (SqlStore): the store object is connected to the backend database and allows for data to be extracted
            and stored.
        defect_key (str): The key of the defect material.
        depth (int): The depth of the defect position in the material.
        expansion_factor (float): The expansion factor of the post MD material compared to the pre MD material.
        defect_index (int): The index of the defect.
        lattice_constant (float): The lattice constant.
    """
    params = MDAbadParameters(
        key=defect_key,
        depth=depth,
        expansion_factor=expansion_factor,
        defect_index=int(defect_index),
        lattice_constant=float(lattice_constant)

    )
    store.save(params)

def CommitAndClose(backend):
    """
    Function to commit changes and close the backend database.
    Args:
        backend : Connection to the database.
    """
    backend.commit()
    backend.close()
