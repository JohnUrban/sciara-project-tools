#!/usr/bin/env python2.7
import sys
import argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""
dnaInSalGld.py is a simple utility for me when thinking about DNA content in larval sal glands.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fly', '-f',
                   type=str, default='sciara',
                   help='''Which fly? Options: sciara, drosophila. Default: sciara.''')

parser.add_argument('--UR', '-u',
                   type=float, default=0.0,
                   help='''Proportion of genome that may be under-replicated (0-1). Default: 0.0.''')

parser.add_argument('--target', '-t',
                   type=float, default=False,
                   help='''Target amount of DNA desired. Default: False. If provided, it will estimate number of pairs needed.''')

parser.add_argument('--ngenomes', '-n',
                   type=float, default=False,
                   help='''Also report the amount of DNA from N genomes. Default: False.''')

parser.add_argument('--nascent', '-N',
                   type=str, default=False,
                   help='''Also report the amount of nascent DNA of size L from N genomes given length of cell cycle (C) and minimum number of origins (O). Default: False.
                        Use --ngenomes/-n to specify the number of cells/genomes. Specify other params here by: --nascent L,C,O.
                        C can be specified in hours or minutes by putting h or m at the suffix of C.
                        This assumes forks move at 1.5 kb/minute.''')

args = parser.parse_args()


fullyRep = 1.0 - args.UR

bpWeight = 1.07935e-12
ntWeight = 0.5*bpWeight
print "Nt weight", ntWeight, "ng"
print "Bp weight", bpWeight, "ng"

if args.fly == 'sciara':
    ngenomesPerCell = 2**13
    nCellsPerPair = 512
    genomeSize = 292e6
elif args.fly == 'drosophila':
    ngenomesPerCell = 2**10
    nCellsPerPair = 200
    genomeSize = 150e6

ngenomesPerPair = nCellsPerPair*ngenomesPerCell

genomeWeight = genomeSize * bpWeight

dnaWeightPerCell = genomeWeight * ngenomesPerCell * fullyRep ## UR factored in at this step

dnaWeightPerPair = dnaWeightPerCell * nCellsPerPair


print "Single Genome:", genomeWeight, "ng"
print "In each cell:", dnaWeightPerCell, "ng"
print "In each pair:", dnaWeightPerPair, "ng"
print "In 10 pairs:", 10*dnaWeightPerPair, "ng"
print "In 20 pairs:", 20*dnaWeightPerPair, "ng"
print "In 20 pairs given 1% yield:", 20*dnaWeightPerPair*0.01, "ng"
print "In 20 pairs given 2% yield:", 20*dnaWeightPerPair*0.02, "ng"
print "In 20 pairs given 5% yield:", 20*dnaWeightPerPair*0.05, "ng"
print "In 20 pairs given 10% yield:", 20*dnaWeightPerPair*0.1, "ng"
print "In 20 pairs given 25% yield:", 20*dnaWeightPerPair*0.25, "ng"
print "In 20 pairs given 50% yield:", 20*dnaWeightPerPair*0.5, "ng"

if args.target:
    print
    print "Minimum number of pairs needed assuming stated percent yields:"
    ans = round(args.target/(dnaWeightPerPair*1.0),1)
    print "100% Yield:", ans, "pairs or", ans/20, "batches of 20."
    ans = round(args.target/(dnaWeightPerPair*0.5),1)
    print "50% Yield:", ans, "pairs or", ans/20, "batches of 20."
    ans = round(args.target/(dnaWeightPerPair*0.25),1)
    print "25% Yield:", ans, "pairs or", ans/20, "batches of 20."
    ans = round(args.target/(dnaWeightPerPair*0.1),1)
    print "10% Yield:", ans, "pairs or", ans/20, "batches of 20."
    ans = round(args.target/(dnaWeightPerPair*0.05),1)
    print "5% Yield:", ans, "pairs or", ans/20, "batches of 20."
    ans = round(args.target/(dnaWeightPerPair*0.01),1)
    print "1% Yield:", ans, "pairs or", ans/20, "batches of 20."

if args.ngenomes:
    print
    print "There is", genomeWeight*args.ngenomes, "ng DNA from", args.ngenomes, "genomes..."
    print "It would take a minimum of", args.ngenomes/ngenomesPerPair, "sg pairs to get this many genomes."

    ## NEED TO CHECK MATH BELOW
    if args.nascent:
        ntWeight = 0.5 * bpWeight
        params = args.nascent.split(',')
        ns_len = int(params[0])
        if params[1].endswith('h'):
            C = float(params[1][:-1])*60
        elif params[1].endswith('m'):
            C = float(params[1][:-1])
        else:
            C = float(params[1])*60
        num_oris = int(params[2])
        fkspd = 1500.0
        num_ns_per_ori = 2 ## could be 4 if you want
        
        ns_life = ns_len/fkspd
        ns_prop_of_cc = ns_life / C

        print "Amount of NS:", args.ngenomes * num_oris * num_ns_per_ori * ns_len * ntWeight * ns_prop_of_cc, "ng"
        print "Amount of NS - max but impossible:", args.ngenomes * num_oris * num_ns_per_ori * ns_len * ntWeight, "ng"


	n_forks = num_oris * 2
	print
	print "Assuming 100% S-phase cells (multiply by proportion in S-phase to get more realistic numbers):"
	print "Okazaki fragment weight (assuming all forks going) = n_forks_per_genome * n_OF_per_fork * nt_weight * n_genome"
	print "Using", args.ngenomes, "genomes (i.e.", 0.5*args.ngenomes, "diploid cells)."
	for i in [0.1, 0.25, 0.5, 0.75,1,2,3,4]:
		print "Max amount of Okazaki fragments assuming", i, "per fork:", n_forks * i * ntWeight * args.ngenomes, "ng"

	
