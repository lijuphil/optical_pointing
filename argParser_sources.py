from __future__ import division, print_function

import os
import argparse

class customArgParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()

def printArgs(args):
    print( 'opFilename: {}'.format(args.pltDir))
    if args.verbose:
        print('verbosity level: {0}'.format(args.verbose))

def get_arguments():
    ### Command line argument parser ###
    parser = customArgParser(description='ArgParser for pointing_sources.py.',fromfile_prefix_chars='@',add_help=False)#,allow_abbrev=False)
    
    parser.add_argument('-sl', '--source_list', type=str, default=None,
                        help='Path to the file listing the sources.')
    parser.add_argument('-pl', '--platform', type=str, default=None,
                        help='Platform, pathfiner or testbed')
    parser.add_argument('-el', '--elevation', nargs='+', type=float, default=[30.0,65.0],
                        help='[min elevation, max elevation] in degrees. Default = [30.0, 65.0]')
    parser.add_argument('-az', '--azimuth_span', type=float, default=18.0,
                        help='Azimuth width in degrees. Default azimuth_span = 18.0 degrees/')
    parser.add_argument('-p', '--plot', nargs='?', default=False,
                        help='Enter True for plots.')

    parse_args, _ = parser.parse_known_args()
    
    # manually adding help argument at the last second ...
    parser.add_argument('-h','--help',action='help',default='==SUPPRESS==',
                        help='show this help message and exit')
    # ... and parse command line now
    parser.parse_args(namespace=parse_args)

    return parse_args
