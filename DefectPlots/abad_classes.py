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

class AbadParameters(httk.HttkObject):

    @httk.httk_typed_init({'key': int, 'depth': float, 'expansion_factor': float, 'defect_index': int})
    def __init__(self, key, depth, expansion_factor, defect_index):
        self.key = key
        self.depth = depth
        self.expansion_factor = expansion_factor
        self.defect_index = defect_index


class Delta(httk.HttkObject):

    @httk.httk_typed_init({'host': str, 'dopant': str, 'defect': str, 'key': int, 'delta': float})
    def __init__(self, host, dopant, defect, key, delta):
        self.host = host
        self.dopant = dopant
        self.defect = defect
        self.key = key
        self.delta = delta