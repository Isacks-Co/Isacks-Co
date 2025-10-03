import argparse
import sys
import ast

import argparse

class SingleMetavarHelpFormatter(argparse.HelpFormatter):
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
        self.argparser.add_argument("-TS", "--Timestep", metavar="<TIMESTEP>", type=float, help="Timestep as a float")
        self.argparser.add_argument("-µ", "--Friction", metavar="<FRICTION>", type=float, help="Friction coefficent as a float (For NVT)")
        self.argparser.add_argument("-C","--Compressibility", metavar="<COMPRESSIBILITY>", type=float, help="Compressibility as a float (NPT)")
        self.argparser.add_argument("-SI", "--Sample_interval", metavar="<INTERVAL>", type=int, help="Sample frequency as an integer" )
        self.argparser.add_argument("-S", "--Supercells", metavar="<SUPERCELL>", type=self.parseList, help="Repetition of input cell e.g [3,3,3], use [1,1,1] for only unit cell")
        self.argparser.add_argument("-O", "--Output_file", metavar="<PATH>", type=str, help="Path to where the output file will be written")
        self.argparser.add_argument("-N", "--Number_of_steps", metavar="<NUMBER_OF_STEPS>", type=str, help="Total number of timesteps as an integer")
        

    def parseList(self, arg):
        """Help function to make sure list arguments are working."""
        try:
            value = ast.literal_eval(arg)
            if type(value) == list:
                return value
            else:
                  raise argparse.ArgumentTypeError("Input must be a list")
        except Exception as e:
            raise argparse.ArgumentTypeError(f"Invalid list: {e}")
    
if __name__ == "__main__":
    parser = InputParser()
    args = parser.argparser.parse_args()
    settings_dict = vars(args)
    print(settings_dict)
