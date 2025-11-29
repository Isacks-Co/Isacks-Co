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

from httk.atomistic import Structure, StructureTag
from httk.core.vectors import FracVector


class newDefectInfos(httk.HttkObject):
    #Class to primary search for info on the specific defect, contains info on host- and defect names, types etc
    #Might not be as useful to us, since we probably only need newDefectCell

    @httk.httk_typed_init({'key': int, 'host_name': str, 'defect_name': str, 'defect_type': str, 'configuration': str, 'vacancy': bool,
                           'substitutional': bool, 'interstitial': bool}, index=['key', 'defect_name',  'defect_type', 'configuration',
                                                                                 'vacancy', 'substitutional', 'interstitial'])
    def __init__(self, key, host_name, defect_name,  defect_type, configuration, vacancy, substitutional, interstitial):
        #I removed some of the parameters, they seem redundant, check sqlite if you want them back
        self.key = key
        self.host_name = host_name
        # unique name for the defect
        self.defect_name = defect_name
        # and where in the supercell it is localated
        self.defect_type = defect_type
        # unique hash for the defect
        self.configuration = configuration
        # defect info
        self.vacancy = vacancy
        self.substitutional = substitutional
        self.interstitial = interstitial

class newDefectCell(httk.HttkObject):
    #Contains the structural information on the defect
    @httk.httk_typed_init({'host_name': str, 'defect_structure': Structure, 'defect_name': str, 'defect_types': [str],
                           'defect_positions': (FracVector, 0, 3), 'key': int },
                          index=['defect_name', 'key'])

    def __init__(self, host_name, defect_structure, defect_name, defect_types, defect_positions, key):
        self.host_name = host_name
        self.defect_structure = defect_structure
        # unique name for the defect
        self.defect_name = defect_name
        # what defect is it
        self.defect_types = defect_types
        # and where in the supercell it is localated
        self.defect_positions = defect_positions
        # unique hash for the defect
        self.key = key
        # priority for the defect
        # general description of the defect

class newScreenResult(httk.Result):

    @httk.httk_typed_init({'defect_key': int, 'computation': httk.Computation,   'defect_folder_name': str, 'total_energy_coarse': float, 'smallest_abs_estimate': float, 'smallest_abs_transition_intensity': float, 'smallest_abs_partial_density_difference': float, 'smallest_abs_polarization_x': float, 'smallest_abs_polarization_y': float, 'smallest_abs_polarization_z': float, 'excitation': str, 'abs_estimate': float, 'abs_transition_intensity': float, 'abs_partial_density_difference': float, 'abs_polarization_x': float, 'abs_polarization_y': float, 'abs_polarization_z': float, 'abs_radiative_lifetime': float, 'ZPL_converged': bool, 'ZPL_estimate': float, 'ZPL_transition_intensity': float, 'ZPL_partial_density_difference': float, 'ZPL_polarization_x': float, 'ZPL_polarization_y': float, 'ZPL_polarization_z': float, 'ZPL_radiative_lifetime': float, 'max_relaxation': float, 'average_relaxation': float, 'delta_R': float, 'delta_Q': float, 'huang_rhys': float, 'debye_waller': float, 'omega': float, 'time': float}, index=['defect_key'])
    def __init__(self, defect_key, computation, defect_folder_name, total_energy_coarse, smallest_abs_estimate, smallest_abs_transition_intensity, smallest_abs_partial_density_difference, smallest_abs_polarization_x, smallest_abs_polarization_y, smallest_abs_polarization_z, excitation, abs_estimate, abs_transition_intensity, abs_partial_density_difference, abs_polarization_x, abs_polarization_y, abs_polarization_z, abs_radiative_lifetime, ZPL_converged, ZPL_estimate, ZPL_transition_intensity, ZPL_partial_density_difference, ZPL_polarization_x, ZPL_polarization_y, ZPL_polarization_z, ZPL_radiative_lifetime, max_relaxation, average_relaxation, delta_R, delta_Q, huang_rhys, debye_waller, omega, time):
        self.defect_key = defect_key
        self.computation = computation
        self.defect_folder_name = defect_folder_name
        self.total_energy_coarse = total_energy_coarse

        # relaxation from ground to excited geometry
        self.max_relaxation = max_relaxation
        self.average_relaxation = average_relaxation






class newDelta(httk.HttkObject):
    """

    """
    @httk.httk_typed_init({'host': str, 'dopant': str, 'defect': str, 'key': int, 'delta': float})
    def __init__(self, host, dopant, defect, key, delta):
        self.host = host
        self.dopant = dopant
        self.defect = defect
        self.key = key
        self.delta = delta