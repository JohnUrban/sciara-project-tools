#!/usr/bin/env python2.7
import sys, argparse, pybedtools, scipy.stats
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description="""

 Given gcmedian.txt and signal file,
 Generate updated signal file corrected for gcfemedian...

 NOTE: if gcfmedian obtained on FE can only use on FE signal.
 If gcmedian done on SPMR, can only use on SPMR signal.
 Etc.

 NOTE2: signal bins being corrected need to be same bins GC analyzed in.
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--signal', '-a', '-i', '-f',
                   type= str,
                   help='''Path to signal bedGraph that has cols: chr, start, end, GCprop, signal-to-correct.Can tell algo what cols to look for.''')
parser.add_argument('--scalefactors', '-b', 
                   type=str, 
                   help='''Path to gc median parameter table.''')

parser.add_argument('--gccol', '-g', 
                   type=int, default=6,
                   help='''Column GC proportions found in. Default = 6''')
parser.add_argument('--signalcol', '-s', 
                   type=int, default=4,
                   help='''Column signal found in. Default = 4''')
parser.add_argument('--chrcol', '-c', 
                   type=int, default=1,
                   help='''Column chr name found in. Default = 1''')
parser.add_argument('--startcol', '-S', 
                   type=int, default=2,
                   help='''Column start coordinate found in. Default = 2''')
parser.add_argument('--endcol', '-E', 
                   type=int, default=3,
                   help='''Column end coordinate found in. Default = 3''')
parser.add_argument('--mean', '-M', 
                   action='store_true', default=False,
                   help='''Subtract GC mean from signal instead of GC median.''')
args = parser.parse_args()


gccol = args.gccol-1
sigcol = args.signalcol-1
chrcol = args.chrcol-1
startcol = args.startcol-1
endcol = args.endcol-1

gcsub = 0
if args.mean:
    gcsub = 1

statdict = defaultdict(list)
with open(args.scalefactors) as table:
    for row in table:
        row = row.strip().split()
        statdict[int(row[0])] = [float(e) for e in row[1:4]]

with open(args.signal) as table:
    for row in table:
        row = row.strip().split()
        gc = int(100.0*float(row[gccol]))
        sig = float(row[sigcol])
        newsig = sig - statdict[gc][gcsub]
        out = [row[chrcol], row[startcol], row[endcol], newsig]
        print ("\t").join([str(e) for e in out])
        


## could also use NS-seq medians to correct its own signal...
## could do that all in the first step -- either to treatment signal, FE signal, or both...
## if do raw signa or raw spmr -- will still need to normalize for copy/number -- and/or do local lambda peak calling
