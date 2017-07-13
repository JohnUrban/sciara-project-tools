#!/usr/bin/env python2.7

import argparse, os, sys

parser = argparse.ArgumentParser(description="""

    Takes in paths to files -- can even do wildcards -- returns abs paths for each.
    
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('files', metavar='FILES', nargs='+',
                               help='The input files - as many as you want.')

args = parser.parse_args()


if len(args.files) == 1 and (args.files[0] == '-' or args.files[0] == 'stdin'):
    args.files = []
    for fpath in sys.stdin:
        args.files += fpath.split()

for fpath in args.files:
    print os.path.abspath(fpath)
