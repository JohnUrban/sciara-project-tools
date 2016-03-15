#!/usr/bin/env python
import sys
import numpy as np
import pandas
import matplotlib.pyplot as plt
import argparse
import re
from collections import defaultdict
import pysam


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in list, computes asm stats including N50 and NG50 values.

    E = sum(contig_size^2)/Genome size. E-size is from the GAGE paper (Salzberg et al,2012, Genome Research). The E-size is designed to answer the question: If you choose a location (a base) in the reference genome at random, what is the expected size of the contig or scaffold containing that location?

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group()
parser_input.add_argument('-i', "--inputfile",
                   type= str, default="-",
                   help='''Input file.''')
parser_input.add_argument('--cmdline', '-c',
                   type= str, default=False,
                   help='''Input list of numbers on cmd line (comma-separated). ''')

parser.add_argument('-k', "--colnum",
                   type=int, default=2,
                   help='''Column number (1-based) to compute n50 on from Input file. Default is first column.''')

parser.add_argument('-x', "--x",
                   type=str, default="50",
                   help='''Give comma-separated X values for NX function -- i.e. 50 for N50. Default=25,50,75''')

parser.add_argument('-G', "--genomesize",
                   type=str, default="210000000,350000000",
                   help='''Produce NG statistics and (if applicable) E-size with some specified genome size (default is G = assembly size = sum(contigs)). Supply comma-separated integer values for genome sizes to try.''')


args = parser.parse_args()


##############################################################################
''' FUNCTIONS '''
##############################################################################

def NX(l, x=[25,50,75], G=False):
        """
        Returns NX for all x for a list of numbers l.
        Default: N25, N50, N75
        Assumes all values in list x are between 0 and 100.
        Interpretation: When NX = NX_value, X% of data (in bp) is contained in reads at least NX_value bp long.
        """
        ## assumes both l and x are sorted
	if isinstance(l, list) and isinstance(x, list) and G:
            l = l[:]
            x = x[:]
            nxsum = 0
            nxvalues = {e:0 for e in x}
            for e in x:
                    xpct = G*e/100.0
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


##############################################################################
''' PROCESS ARGS '''
##############################################################################


## is input from file format?
if args.inputfile:
    # is file format from stdin or from specified location?
    if args.inputfile == "stdin" or args.inputfile == "-":
        connection = sys.stdin
    else:
        connection = open(args.inputfile, 'r')
        
    # require the file to have at least 1 column to work on
    assert args.colnum > 0
    
    ## attempt to extract list of sizes from input
    l = []
    try:
        for line in connection:
            line = line.strip().split()
            l.append(float(line[args.colnum-1]))
    except IndexError:
        print "This file does not have this many columns... please provide tabdelim file"
        quit()

## or is input command-line list format?       
elif args.cmdline:
    # if list given at command-line, process it
    l = [float(e) for e in args.cmdline.split(",")]



## sort sizes
l.sort()

## get X values for NX stats and sort
x = [float(e) for e in args.x.split(",")]
x.sort()

## get genome size values and sort
G = [float(e) for e in args.genomesize.split(",")]
G.sort()

## Get N contigs
print "Number contigs:", len(l)

## Get assembly size
A = sum(l)
print "Assembly size:", A

## Get max contig size:
MAX = max(l)
print "Max contig size:", MAX

## Get min contig size:
MIN = min(l)
print "Min contig size:", MIN

## Get mean contig size
MEAN = np.mean(l)
print "Mean contig size:", MEAN

## Get median contig size
MEDIAN = np.median(l)
print "Median contig size:", MEDIAN

## Get NX values
nxvalues = NX(l,x,G=A)
for e in x:
    print "Contig N%s\t%d" % (str(e), nxvalues[e])


## expected value given assembly size
E = e_size(l,G=A)
print "E size (G=%d) = %d" % (A, E)

## Get NGX values
for g in G:
    nxvalues = NX(l,x,G=g)
    for e in x:
        print "Contig N%s (G=%d)\t%d" % (str(e), g, nxvalues[e])

## get expected sizes given genome size values
for g in G:
    E = e_size(l,G=g)
    print "E size (G=%d) = %d" % (g, E)





