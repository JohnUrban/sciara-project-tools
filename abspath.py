#!/usr/bin/env python2.7

import argparse, os, sys

parser = argparse.ArgumentParser(description="""

    Takes in paths to files -- can even do wildcards -- returns abs paths for each.
    
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('files', metavar='FILES', nargs='+', help='The input files - as many as you want.')
parser.add_argument('-b', '--basename', action='store_true', default=False, help='''Return basename instead.''')
parser.add_argument('-s', '--split', action='store_true', default=False, help='''Return path to dir and basename separately/split.''')

args = parser.parse_args()


if len(args.files) == 1 and (args.files[0] == '-' or args.files[0] == 'stdin'):
    args.files = []
    for fpath in sys.stdin:
        args.files += fpath.split()

for fpath in args.files:
    if args.basename:
        print os.path.basename(fpath)
    elif args.split:
        print ("\t").join( os.path.split(os.path.abspath(fpath)) )
    else: #abspath
        print os.path.abspath(fpath)
