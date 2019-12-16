#!/usr/bin/env python2.7
import argparse, sys
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
Description: Takes in Fasta.
Checks beginning and end of each sequence to ensure it does not start with N (or w/e specified).
''')
parser.add_argument("fasta")
parser.add_argument("-G", "--gapchar", type=str, default='N', help='''The single letter gap character. Default =  N.''')
args = parser.parse_args()


# Open FASTA, search for masked regions, print in BED3 format
N = args.gapchar + '+' ## 'N+'
handle = sys.stdin if args.fasta in ('-','stdin') else open(args.fasta)
pattern = N

## HEADER
contents = ['#', 'name', 'length', 'leading_N', 'trailing_N', 'first_N_pos', 'last_N_pos']
print '\t'.join(contents)

## RECORDS
for record in SeqIO.parse(handle, "fasta"):
    length = len(record.seq)
    allmatch = [e for e in re.finditer(pattern, str(record.seq))]
    if len(allmatch) > 0:
        first = allmatch[0].start()
        last = allmatch[-1].end()
        start = 'yes' if first == 0 else 'no'
        end = 'yes'if last == length else 'no'
        contents = [record.id, length, start, end, first, last]
    else:
        contents = [record.id, length, 'no','no', 'na', 'na']
    print '\t'.join([str(e) for e in contents])

    
if args.fasta not in ('-','stdin'):
    handle.close()
