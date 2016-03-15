#!/usr/bin/env python
import sys, argparse
import numpy as np
from Bio import SeqIO
from helper_functions import *
##from joblib import Parallel, delayed


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in fasta or fastq assembly file. 
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument("-fq", "--fastq",
                   type= str, default=False,
                   help='''Input fastq file.''')
parser_input.add_argument("-fa", "--fasta",
                   type= str, default=False,
                   help='''Input fasta file.''')
##parser.add_argument("-t", type=int, default=1,
##                    help=''' Number of threads. Default: 1.''')
args = parser.parse_args()


contig_qvs = []
contig_lengths = []
allqv = []
if args.fasta:
    args.fastx = args.fasta
    fastx = "fasta"
elif args.fastq:
    args.fastx = args.fastq
    fastx = "fastq"


def get_info(contig):
    length = len(contig)
    n_lower = case_counter(str(contig.seq))
    n_upper = length - n_lower
    pct_lower = 100.0*n_lower/length
    pct_upper = 100.0*n_upper/length
    sys.stdout.write(("\t").join([contig.name, str(length), str(n_lower), str(n_upper), str(pct_lower), str(pct_upper)])+"\n")


for contig in SeqIO.parse(args.fastx, "fastq"):
    get_info(contig)
##    length = len(contig)
##    n_lower = case_counter(str(contig.seq))
##    n_upper = length - n_lower
##    pct_lower = 100.0*n_lower/length
##    pct_upper = 100.0*n_upper/length
##    sys.stdout.write(("\t").join([contig.name, str(length), str(n_lower), str(n_upper), str(pct_lower), str(pct_upper)])+"\n")



##Parallel(n_jobs=args.t)(delayed(get_info)(contig) for contig in SeqIO.parse(args.fastx, "fastq"))

