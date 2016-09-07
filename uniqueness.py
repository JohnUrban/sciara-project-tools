#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict

## DEFINE PARSER ARGUMENTS
parser = argparse.ArgumentParser(description="""

USAGE RECOMMENDATION:
    uniqueness.py file.fasta > uniqueness.wig
    or
    uniqueness.py --bdg file.fasta > uniqueness.bedGraph
    or
    uniqueness.py --bdg file.fasta | sortBed -i - > uniqueness.bedGraph
    
    (Change kmer size as necessary with -k)
    
DESCRIPTION
    
    This is slow and needs a lot of memory. Don't use it... unless you need to...
    
    Take in fasta.
    Count number of times each kmer occurs.
    Assign each position containing kmer K a uniques scor of 1/|K|
    ..where |K| is the number of times kmer K appeared.


    TODO:
    return in wig format option.
    In wig, for each sequence/chromosome/contig, one specifies a header followed by 1 score per line -- e.g.:
    fixedStep chrom=000000F|quiver|pilon	start=0	step=10000 span=20000
    1
    2
    3
    4
    5
    6
    5
    4
    3
    2
    1
    ...
    
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('fasta', metavar='fasta', nargs='+',
                   type= str, 
                   help='''Path(s) to fasta file(s).
                        Can handle more than one, though it might not be recommended to use more than one.
                        If provide more than 1, then it still assumes all sequence names are uniqe.''')
parser.add_argument('-k', '--kmersize', type=int, default=50,
                    help='''Default = 50.''')

parser.add_argument('-b', '--bdg', action='store_true',
                    help='''Default output is wig. This outputs bedGraph''')

parser.add_argument('-v', '--verbose', action='store_true')

##parser.add_argument('--step', type=str, default='1',
##                    help='''step size for fixed step wigggle format. default = 1. DO NOT CHANGE THIS.''')
##parser.add_argument('--span', type=str, default='1',
##                    help='''span for fixed step wigggle format. default = 1. DO NOT CHANGE THIS.''')

args = parser.parse_args()


##kmerlocs = defaultdict(list) ## list of tuples...
kmercounts = defaultdict(int)
step = '1'
span = '1'
k=args.kmersize

def bedgraph(args, kmercounts,k):
    for f in args.fasta:
        for fa in SeqIO.parse(f, "fasta"):
            seq = str(fa.seq).upper()
            name = str(fa.name)
            for i in range(len(seq)-k+1):
                u = 1.0/kmercounts[seq[i:i+k]]
                print ("\t").join([str(e) for e in [name, i, i+1, u]])

def wiggle(args, kmercounts,k):
    for f in args.fasta:
        for fa in SeqIO.parse(f, "fasta"):
            seq = str(fa.seq).upper()
            name = str(fa.name)
            print "fixedStep chrom=" + name + " start=0 step=" + step + " span=" + span
            for i in range(len(seq)-k+1):
                u = 1.0/kmercounts[seq[i:i+k]]
                print str(u)


if args.verbose:
    sys.stderr.write("Counting kmers....\n")
for f in args.fasta:
    if args.verbose:
        sys.stderr.write("Counting kmers...."+f+"\n")
    for fa in SeqIO.parse(f, "fasta"):
        seq = str(fa.seq).upper()
        if args.verbose:
            sys.stderr.write("Counting kmers...."+f+"...Working on "+str(fa.name)+"...\n")#Len = "+str(len(seq))+" ...num unique kmers seen before this step = "+str(len(kmercounts.keys()))+"....\n")
        for i in range(len(seq)-k+1):
            kmercounts[seq[i:i+k]] += 1

if args.verbose:
    sys.stderr.write("Returning output....\n")
if args.bdg:           
    bedgraph(args,kmercounts,k)
else: ##wiggle
    wiggle(args,kmercounts,k)
