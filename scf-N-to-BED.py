#!/usr/bin/env python2.7
import argparse, sys
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
Description: Takes in Fasta.
Outputs BED of N-gap locations.
Use '-' or 'stdin' if coming from stdin.
Use -l/--length to control minimum gap length (default 25).
Use -r/--recognition to allow regex to expand gaps through recognition sequence islands dispersed in N-gaps (e.g. from optical maps).
When using -r provide only the recognition sequence. For now it is case-sensitive.
By default N-gaps are simply where Ns start to the first non-N.
''')
parser.add_argument("fasta")
parser.add_argument("-l", "--length", type=int, default=25, help='''Do not report as gap if less than this length. Default=25.''')
parser.add_argument("-r", "--recognition", type=str, default=False, help='''Allow N-gaps to encompass occurences of given recognition sequence. Default=False.''')
parser.add_argument("-G", "--gapchar", type=str, default='N', help='''The single letter gap character. Default =  N.''')
args = parser.parse_args()


# Open FASTA, search for masked regions, print in BED3 format
N = args.gapchar + '+' ## 'N+'
handle = sys.stdin if args.fasta in ('-','stdin') else open(args.fasta)
pattern = N
if args.recognition:
    pattern = N + '((' + args.recognition + '){0,1}' + N + ')*' ## 1+ Ns and optionally 0 or more (recseq, Ns) 
for record in SeqIO.parse(handle, "fasta"):
    for match in re.finditer(pattern, str(record.seq)):
        if match.end()-match.start()  >= args.length:
            print ("\t").join([str(e) for e in (record.id, match.start(), match.end())])

if args.fasta not in ('-','stdin'):
    handle.close()
