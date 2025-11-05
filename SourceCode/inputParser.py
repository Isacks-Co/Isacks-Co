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

from logger import loggerSetup

import argparse
import sys
import ast
import logging

class SingleMetavarHelpFormatter(argparse.HelpFormatter):
    """Formatter to only get one metavar"""
    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)

        opts = ', '.join(action.option_strings)
        metavar = action.metavar or (action.type.__name__.upper() if action.type else None)

        if metavar:
            return f'{opts} {metavar}'
        else:
            return opts

class InputParser():
    def __init__(self, args=sys.argv):
        self.argparser = argparse.ArgumentParser(prog="MolecularDynamics.py", 
                                                 usage='%(prog)s <Atomic Configuration> <settings> [options]',
                                                 description="Molecular dynamics program developed by Isacks & Co", 
                                                 formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, max_help_position=50, width=100))
        self.addArguments()
        self.input_args = args
        self.args = vars(self.argparser.parse_args(args[1:]))
        
        loggerSetup(self.args["debug"])

    def addArguments(self):
        """ Whenever a new setting is added, add the corresponding flag here."""
        #Positional arguments
        self.argparser.add_argument("input_structure", help="Path for atomic configuration")
        self.argparser.add_argument("input_settings", help="Path to settings file")
        #Options (flags), to add a new flag make sure to follow the format below with -<abbrev.>, --<Full name>, etc.
        self.argparser.add_argument("-E", "--Ensemble", metavar="<ENSEMBLE>", type=str, help="Ensemble (NVE, NPT, NVT)")
        self.argparser.add_argument("-T","--Temperature", metavar="<TEMPERATURE>", type=float, help="Temperature in K")
        self.argparser.add_argument("-P","--Pressure", metavar="<PRESSURE>", type=float, help="Pressure in Pa")
        self.argparser.add_argument("-POT","--Potential" , metavar="<POTENTIAL>", type=str, help="Potential as a string (EMT, LJ, MACE)")
        self.argparser.add_argument("-TS", "--Timestep", metavar="<TIMESTEP>", type=float, help="Timestep as a float (fs)")
        self.argparser.add_argument("-µ", "--Friction", metavar="<FRICTION>", type=float, help="Friction coefficent as a float (For NVT)")
        self.argparser.add_argument("-TD", "--Tdamp", metavar="<TDAMP>", type=float, help="Tdamp as a float (For NPT)")
        self.argparser.add_argument("-PD", "--Pdamp", metavar="<PDAMP>", type=float, help="Pdamp as a float (For NPT)")
        self.argparser.add_argument("-SI", "--Sample_interval", metavar="<INTERVAL>", type=int, help="Sample frequency as an integer" )
        self.argparser.add_argument("-S", "--Supercells", metavar="<SUPERCELL>", type=self.parseList, help="Repetition of input cell e.g [3,3,3], use [1,1,1] for only unit cell")
        self.argparser.add_argument("-O", "--Output_file", metavar="<PATH>", type=str, help="Path to where the output file will be written")
        self.argparser.add_argument("-N", "--Number_of_steps", metavar="<NUMBER_OF_STEPS>", type=str, help="Total number of timesteps as an integer")
        self.argparser.add_argument("--debug", action="store_true", help="Debug")
        

    def parseList(self, arg):
        """Help function to make sure list arguments are working."""
        try:
            value = ast.literal_eval(arg)
            if type(value) == list:
                return value
            else:
                  raise argparse.ArgumentTypeError("Input must be a list")
        except Exception:
            raise argparse.ArgumentTypeError(f"Invalid list: {arg}")
if __name__ == "__main__":
    parser = InputParser()
    args = parser.argparser.parse_args()
    settings_dict = vars(args)
    print(settings_dict)
