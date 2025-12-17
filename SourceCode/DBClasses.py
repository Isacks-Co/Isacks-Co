###!/usr/bin/env python
#
#    Automatic Defect Analysis and Qualification (ADAQ)
#    Copyright (C) 2016-2021 Joel Davidsson
#    Implemented using the high-throughput toolkit (httk) (see README.md)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function, division

try:
    import httk
except Exception:
    import sys, os.path, inspect

    _realpath = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
    sys.path.insert(1, os.path.join(_realpath, '../..'))
    import httk

from httk.atomistic import Structure


class MDDefectInfos(httk.HttkObject):
    """
    Metadata container for defect classification and search.

    Stores human- and workflow-friendly defect metadata such as defect type,
    configuration identifier, and boolean flags (vacancy/substitutional/interstitial).
    Intended primarily for searching and filtering defects in the database.

    Attributes
    ----------
    key : int
        Unique identifier for the defect entry.
    host_name : str
        Name/identifier of the host material.
    defect_name : str
        Name/identifier of the defect (often includes dopant and site/type).
    defect_type : str
        High-level defect category label.
    configuration : str
        Configuration label (e.g., site index or symmetry-derived tag).
    vacancy : bool
        True if the defect is a vacancy type.
    substitutional : bool
        True if the defect is substitutional.
    interstitial : bool
        True if the defect is interstitial.
    """
    @httk.httk_typed_init(
        {'key': int, 'host_name': str, 'defect_name': str, 'defect_type': str, 'configuration': str, 'vacancy': bool,
         'substitutional': bool, 'interstitial': bool}, index=['key', 'defect_name', 'defect_type', 'configuration',
                                                               'vacancy', 'substitutional', 'interstitial'])
    def __init__(self, key, host_name, defect_name, defect_type, configuration, vacancy, substitutional, interstitial):
        self.key = key
        self.host_name = host_name
        self.defect_name = defect_name
        self.defect_type = defect_type
        self.configuration = configuration
        self.vacancy = vacancy
        self.substitutional = substitutional
        self.interstitial = interstitial


class MDDefectCell(httk.HttkObject):
    """
    Database object for a relaxed defect structure.

    Stores the relaxed defect cell structure associated with a host material
    and a defect name/key.

    Attributes
    ----------
    host_name : str
        Name/identifier of the host material (often includes method tags like "_PBE").
    defect_structure : httk.atomistic.Structure
        Relaxed structure stored in HTTK format.
    defect_name : str
        Defect identifier (often includes dopant and configuration label).
    key : int
        Unique defect key used to join against other tables/results.
    """
    @httk.httk_typed_init({'host_name': str, 'defect_structure': Structure, 'defect_name': str,
                           'key': int},
                          index=['defect_name', 'key'])
    def __init__(self, host_name, defect_structure, defect_name, key):
        self.host_name = host_name
        self.defect_structure = defect_structure
        self.defect_name = defect_name
        self.key = key


class MDScreenResult(httk.Result):
    """
    Store stability delta between interstitial and adatom configurations.

    The delta is defined as:

        delta = E_interstitial_min - E_adatom_min

    Interpretation
    --------------
    - delta < 0: interstitial is more stable
    - delta > 0: adatom is more stable

    Attributes
    ----------
    host : str
        Host material identifier.
    dopant : str
        Dopant identifier.
    defect : str
        Defect name corresponding to the chosen minimum-energy configuration.
    key : int
        Defect key corresponding to `defect`.
    delta : float
        Energy difference in eV (as defined above).
    """
    @httk.httk_typed_init({
        'defect_key': int,
        'computation': httk.Computation,
        'defect_folder_name': str,
        'total_energy_coarse': float,
        'max_relaxation': float,
        'average_relaxation': float
    }, index=['defect_key'])
    def __init__(self, defect_key, total_energy_coarse):
        self.defect_key = defect_key
        self.total_energy_coarse = total_energy_coarse


class MDDelta(httk.HttkObject):
    """
    Store stability delta between interstitial and adatom configurations.

    The delta is defined as:

        delta = E_interstitial_min - E_adatom_min

    Interpretation
    --------------
    - delta < 0: interstitial is more stable
    - delta > 0: adatom is more stable

    Attributes
    ----------
    host : str
        Host material identifier.
    dopant : str
        Dopant identifier.
    defect : str
        Defect name corresponding to the chosen minimum-energy configuration.
    key : int
        Defect key corresponding to `defect`.
    delta : float
        Energy difference in eV (as defined above).
    """
    @httk.httk_typed_init({'host': str, 'dopant': str, 'defect': str, 'key': int, 'delta': float},
                          index = ["key"])
    def __init__(self, host, dopant, defect, key, delta):
        self.host = host
        self.dopant = dopant
        self.defect = defect
        self.key = key
        self.delta = delta

class MDAbadParameters(httk.HttkObject):
    """
    Auxiliary geometric parameters for post-processing and filtering.

    Stores additional values computed from relaxed structures that are used
    in downstream analysis (e.g., defect depth and expansion filtering).

    Attributes
    ----------
    key : int
        Defect key linking these parameters to a defect entry.
    depth : float
        Defect depth metric (typically along the surface normal / z-direction).
    expansion_factor : float
        Relative expansion measure used to detect unstable/"exploded" runs.
    defect_index : int
        Index of the defect atom in the atom list used for depth calculations.
    lattice_constant : float
        Lattice constant metric used for later analysis.

    Notes
    -----
    The ``httk_typed_init`` schema in the decorator should include
    ``lattice_constant`` if it is intended to be stored in the database.
    """
    @httk.httk_typed_init({'key': int, 'depth': float, 'expansion_factor': float, 'defect_index': int},
                          index = ["key"])
    def __init__(self, key, depth, expansion_factor, defect_index, lattice_constant):
        self.key = key
        self.depth = depth
        self.expansion_factor = expansion_factor
        self.defect_index = defect_index
        self.lattice_constant = lattice_constant
