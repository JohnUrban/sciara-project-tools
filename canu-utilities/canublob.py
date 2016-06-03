#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from matplotlib import pyplot as plt
from helper_functions import gc
from plot_functions import scatterplot
from math import log10, sqrt
## DEFINE PARSER ARGUMENTS
parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in maligner files and assembly, reports potential scaffolds of contigs.
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('canu', metavar='canu', nargs='+',
                   type= str, 
                   help='''Path(s) to canu fasta file(s).''')
parser.add_argument('-s', '--saveas', type=str, default=None,
                    help='''Instead of plotting to screen, save as either jpg or pdf by specifying -s name.jpg or -s name.pdf.''')
parser.add_argument('-S', '--scale', type=float, default=1e3,
                    help='''''')
parser.add_argument('-n', '--names', action='store_true', default=False, help='''Plot contig names on points.''')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='''Say stuff whle doing stuff.''')
args = parser.parse_args()



## GRAB INFO
COV = []
GC = []
L = []
names = []
covstats = []
if args.verbose: sys.stderr.write("Starting...\n")
i = 0
j = 0
for canu in args.canu:
    i += 1
    if args.verbose: sys.stderr.write("File %d...\n" % (i))
    for fa in SeqIO.parse(canu, "fasta"):
        j += 1
        if args.verbose and j%10 == 0: sys.stderr.write("Sequence %d..\n" % (j)) 
        desc = fa.description.split()
        nreads = float(desc[2].split("=")[1])
        length = float(desc[1].split("=")[1])
        covstat = float(desc[3].split("=")[1])
        covstats.append(covstat)
        L.append(length)
        cov = nreads/length
        logcov = log10(cov)
        gccontent = gc(str(fa.seq))
        COV.append(logcov)
        GC.append(gccontent)
        names.append(desc[0])


S = [e/args.scale for e in L]

cov = []
gc = []
s = []
for e in sorted(zip(S,COV,GC), reverse=True):
    s.append(e[0])
    cov.append(e[1])
    gc.append(e[2])


pos = float(max(covstats))
pos = pos if pos > 0 else 1
neg = float(min(covstats))
neg = abs(neg) if neg < 0 else 1
covstat2 = []
for e in covstats:
    if e > 0:
        e = e/pos
    elif e < 0:
        e = e/neg
    covstat2.append(e)
covstats = covstat2

print min(covstats), max(covstats)

## PLOT BLOB - NOTE - in future, should plot in order of longest to shortest -- so biggest circles are always on bottom
if not args.names:
    names = False

scatterplot(X=gc, Y=cov, vmin=min(covstats), vmax=max(covstats), size=s, names=names, color=covstats, colorbar=True, alpha=0.7, saveas=args.saveas)#, cmap='seismic') 
       
