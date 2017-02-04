#!/usr/bin/env python
import sys, argparse

parser = argparse.ArgumentParser(description="""
    By John Urban.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-b", "--blastfile",
                   type= str, default=False, required=True,
                   help='''Tab-delim blast output file.''')
parser.add_argument("-k", "--columns",
                   type=str, default="2,13,14,8",
                   help='''1-based, comma-separated column numbers for (in this order) sequence name, start, and end (for the subject), and strand. Default: 2,13,14,8.''')
parser.add_argument("-K", "--othercolumns",
                   type=str, default=False, 
                   help='''1-based Column number for other columns you want to trail.''')

args = parser.parse_args()


cols = [int(e) - 1 for e in args.columns.split(",")]
name = cols[0]
start = cols[1]
end = cols[2]
strand = cols [3]
others = []
if args.othercolumns:
    others = [int(e) - 1 for e in args.othercolumns.split(",")]
if args.blastfile == "-" or args.blastfile =="stdin":
    blastfile = sys.stdin
else:
    blastfile = open(args.blastfile, 'r')

for entry in blastfile:
    entry = entry.strip().split()
    if entry[strand] == "plus":
        print ("\t").join([entry[name], str(int(entry[start])-1), entry[end]] + [entry[e] for e in others])
    elif entry[strand] == "minus":
        print ("\t").join([entry[name], str(int(entry[end])-1), entry[start]] + [entry[e] for e in others])


    
