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
    Save an :class:`DBClasses.MDDefectCell` entry to the database.

    Parameters
    ----------
    store : SqlStore
        Store object connected to the backend database.
    host_name : str
        Name/identifier of the host material.
    defect_name : str
        Name/identifier of the defect configuration.
    defect_key : str or int
        Defect key identifier. Converted to ``int`` before storing.
    structure
        Defect structure object stored in the ``defect_structure`` field.

    Returns
    -------
    None
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
    Save an :class:`DBClasses.MDScreenResult` entry to the database.

    Parameters
    ----------
    store : SqlStore
        Store object connected to the backend database.
    defect_key : str or int
        Defect key identifier. Converted to ``int`` before storing.
    total_energy : float
        Total energy value stored as ``total_energy_coarse``.

    Returns
    -------
    None
    """
    result = MDScreenResult(
        defect_key=int(defect_key),
        total_energy_coarse=total_energy,
    )
    store.save(result)


def saveMDAbadParameters(store, defect_key, depth, expansion_factor, defect_index, lattice_constant):
    """
    Save an :class:`DBClasses.MDAbadParameters` entry to the database.

    Parameters
    ----------
    store : SqlStore
        Store object connected to the backend database.
    defect_key : str or int
        Defect key identifier stored as ``key``.
    depth : float
        Defect depth metric (typically along the z-direction).
    expansion_factor : float
        Relative expansion measure used to detect unstable/"exploded" runs.
    defect_index : int
        Index of the defect atom.
    lattice_constant : float
        Lattice constant metric used in later analysis.

    Returns
    -------
    None
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
    Commit pending changes and close the database backend.

    Parameters
    ----------
    backend
        Database backend connection object providing ``commit()`` and ``close()``.

    Returns
    -------
    None
    """
    backend.commit()
    backend.close()
