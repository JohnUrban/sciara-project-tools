#!/usr/bin/env python2.7
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta")
parser.add_argument("-l", "--length", type=int, default=25, help='''Do not report as gap if less than this length. Default=25.''')
args = parser.parse_args()

# Open FASTA, search for masked regions, print in BED3 format
with open(args.fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('N+', str(record.seq)):
            if match.end()-match.start()  >= args.length:
                print ("\t").join([str(e) for e in (record.id, match.start(), match.end())])
