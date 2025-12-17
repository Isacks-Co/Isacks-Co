#!/usr/bin/env python
#
# Description:
#    This file was created as part of the Automatic Defect Analysis and Qualification (ADAQ) software packages.
#    For additional information about ADAQ, visit: https://defects.anyterial.se/
#
# License:
#    This file is part of the project distributed under the MIT License.
#    Copyright (c) 2025 Joel Davidsson
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

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
