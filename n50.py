#!/usr/bin/env python
import sys
import numpy as np
###import pandas
###import matplotlib.pyplot as plt
import argparse
###import re
from collections import defaultdict
###import pysam


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in list, computes N50 or NX values.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('-i', "--inputfile",
                   type= str, default=False,
                   help='''Input file.''')
parser_input.add_argument('--cmdline', '-c',
                   type= str, default=False,
                   help='''Input list of numbers on cmd line (comma-separated). ''')

parser.add_argument('-k', "--colnum",
                   type=int, default=1,
                   help='''Column number (1-based) to compute n50 on from Input file. Default is first column.''')

parser.add_argument('-x', "--x",
                   type=str, default="25,50,75",
                   help='''Give comma-separated X values for NX function -- i.e. 50 for N50. Default=25,50,75''')

parser.add_argument('-E', "--gageEsize",
                   action="store_true", default=False,
                   help='''E = sum(contig_size^2)/Genome size. E-size is from the GAGE paper (Salzberg et al,2012, Genome Research). The E-size is designed to answer the question: If you choose a location (a base) in the reference genome at random, what is the expected size of the contig or scaffold containing that location? ''')

parser.add_argument('-G', "--genomesize",
                   type=int, default=False,
                   help='''Produce NG statistics and (if applicable) E-size with some specified genome size (default is G=sum(contigs)). Supply integer value for genome size.''')

parser.add_argument('-pctdatagtx',
                   type=str, default=False,
                   help='''Instead of NX values, return pct of data greater than X. Provide X with this flag.''')
parser.add_argument('-pctreadsgtx',
                   type=str, default=False,
                   help='''Instead of NX values, return pct of items (reads, contigs) in list greater than X. Provide X with this flag.''')


args = parser.parse_args()



def NX(l, x=[25,50,75], G=False):
        """
        Returns NX for all x for a list of numbers l.
        Default: N25, N50, N75
        Assumes all values in list x are between 0 and 100.
        Interpretation: When NX = NX_value, X% of data (in bp) is contained in reads at least NX_value bp long.
        """
	if isinstance(l, list) and isinstance(x, list):
		l = sorted(l)
		x = sorted(x)
		if G:
                    total = G
                else:
                    total = sum(l)
                nxsum = 0
                nxvalues = {e:0 for e in x}
		for e in x:
                        xpct = total*e/100.0
                        while nxsum < xpct and l:
                                nxsum += l[-1]
                                lastsize = l.pop()
                        nxvalues[e] = lastsize
                return nxvalues

	else:
		return None

def e_size(l,G=False):
    if G:
        total = G
    else:
        total = sum(l)
    return sum([e**2 for e in l])/float(total)


##def pctDataGtX(l, x):
##        """
##        Returns percent of data on reads in l that is > x for all x
##        """
##        total = sum(l)
##        gtx = 0
##        for e in l:
##                if e > x:
##                        gtx += e
##        return 100.0*gtx/total
##
##def pctReadsGtX(l, x):
##        """
##        Returns percent elements in l that is > x 
##        """
##        total = len(l)
##        gtx = 0
##        for e in l:
##                if e > x:
##                        gtx += 1
##        return 100.0*gtx/total

if args.inputfile:
    if args.inputfile == "stdin" or args.inputfile == "-":
        connection = sys.stdin
    else:
        connection = open(args.inputfile, 'r')
    assert args.colnum > 0
    l = []
    try:
        for line in connection:
            line = line.strip().split()
            l.append(float(line[args.colnum-1]))
    except IndexError:
        print "This file does not have this many columns... please provide tabdelim file"
        quit()
elif args.cmdline:
    l = [float(e) for e in args.cmdline.split(",")]



if args.pctreadsgtx:
        print pctReadsGtX(l, float(args.pctreadsgtx))
elif args.pctdatagtx:
        print pctDataGtX(l, float(args.pctdatagtx))
else: ## NX operation default
        l.sort()
        x = [float(e) for e in args.x.split(",")]
        nxvalues = NX(l,x,G=args.genomesize)
        if args.genomesize:
            out = "NG%s\t%d"
            G = args.genomesize
        else:
            out = "N%s\t%d"
            G = sum(l)
        for e in x:
            print out % (str(e), nxvalues[e])
        if args.gageEsize:
            E = e_size(l,G=args.genomesize)
            print "E size = %d, genome size used = %d" % (E, G)

