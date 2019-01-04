#!/usr/bin/env python2.7
import sys, argparse

parser = argparse.ArgumentParser(description="""
    By John Urban.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-b", "--blastfile",
                   type= str, default=False, required=True,
                   help='''Tab-delim blast output file.''')
parser.add_argument("-k", "--columns",
                   type=str, default="2,13,14,8",
                   help='''1-based, comma-separated column numbers for (in this order) sequence name, start, and end (for the subject), and strand. Default: 2,13,14,8.
NOTE: if strand is NOT part of columns, just provide the first 3 column numbers.''')
parser.add_argument("-K", "--othercolumns",
                   type=str, default=False, 
                   help='''1-based Column number for other columns you want to trail.''')

args = parser.parse_args()


cols = [int(e) - 1 for e in args.columns.split(",")]
name = cols[0]
start = cols[1]
end = cols[2]
try:
    strand = cols [3]
except:
    strand = False
    
others = []
if args.othercolumns:
    others = [int(e) - 1 for e in args.othercolumns.split(",")]
if args.blastfile == "-" or args.blastfile =="stdin":
    blastfile = sys.stdin
else:
    blastfile = open(args.blastfile, 'r')

for entry in blastfile:
    entry = entry.strip().split()
    if (strand != False and entry[strand] == "plus") or (int(entry[start]) < int(entry[end])):
        print ("\t").join([entry[name], str(int(entry[start])-1), entry[end]] + [entry[e] for e in others])
    elif (strand != False and entry[strand] == "minus") or (int(entry[end]) < int(entry[start])):
        print ("\t").join([entry[name], str(int(entry[end])-1), entry[start]] + [entry[e] for e in others])


    
