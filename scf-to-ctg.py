#!/usr/bin/env python2.7
import sys
import argparse
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Take in FASTA file.

Break sequence at Ns.

Return FASTA file with each scf broken into sequence between Ns (gaps).

Naming will be:
$(scfname)_ctg$(C)_len$(L)bp_after$(N)Ns

Where:
    scfname     = original fasta entry name
    C           = contig number along scaffold
    L           = contig length (after breaking off from scaffold)
    N           = number of Ns between this and previous contig

e.g.
scaff1_ctg1_len10153bp_after0Ns

Given n contigs in a scaffold, there are n-1 intervening regions/gaps.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fasta', '-f',
                   type=str, required=True,
                   help='''Path to input fasta file''')

parser.add_argument('--minctglen', '-m',
                   type=int, default=0,
                   help='''After breaking scaffolds on Ns, minimum contig size to report.''')

parser.add_argument('--excludedseqs', '-e',
                   type=str, default='',
                   help='''Comma-separated list of (usually) short (e.g. 6-7 bp) sequences that should be excluded as individual contigs.
Use case might be the restriction sites left inside large gap regions from BioNano scaffolding. These will be counted as Ns.''')

args = parser.parse_args()

excludedseqs = args.excludedseqs.split(",")


for fa in SeqIO.parse(args.fasta, 'fasta'):
    scf = str(fa.seq)
    ctglist = scf.split('N')
    ctgnum = 0
    Ncnt = 0
    for candidate in ctglist:
        if candidate:
            if len(candidate) > args.minctglen and candidate not in excludedseqs:
                ctgnum += 1
                length = len(candidate)
                print ">" + fa.name + "_ctg" + str(ctgnum) + "_len"+str(length) + "bp_after" + str(Ncnt) + "Ns"
                print candidate
                Ncnt = 1 ## If start at 0, it seems to report number of Ns = true_num - 1
            else:
                Ncnt += len(candidate) + 1 ## If start at 0, it seems to report number of Ns = true_num - 1
        else:
            Ncnt += 1
