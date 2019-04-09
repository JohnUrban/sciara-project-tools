#!/usr/bin/env python2.7
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
    Description: Takes in Fasta of multiple sequences.
    Outputs Fasta w/ those sequences concatenated together with N-gaps AND BED coordinates for each sequence on the pseudo-scaffold.
    IF you map things to the pseudo-scaffold and have those locations in BED format,
    THEN you can get the BED coordinates for the input_sequences by:
        $ intersectBed -wo -a output.pseudo-scaffold.bed -b features-mapped-to-pseudo-scaffold.bed | awk 'OFS="\t" {print $4, $8-$2,$9-$2}' $

    ...where output.pseudo-scaffold.bed is the BED output by this script.
        

    This script currently builds the pseudo-scaffold and BED in memory, and writes it all at once at the end.
    For very large files, it could be re-written to append each BED entry and sequence to the growing FASTA and BED files as it goes.

    See also: fastxPaste.py
''')
parser.add_argument("fasta")
parser.add_argument("-l", "--length", type=int, default=10000, help='''Make N-gaps of this length between sequences. Default=10000. The default is large to help prevent spanning long-read alignments.''')
parser.add_argument("-n", "--name", type=str, default='pseudo_scaf', help='''Name of the pseudo-scaffold. Default: pseudo_scaf.''')
parser.add_argument("-d", "--commentdelim", type=str, default=',', help='''For fasta headers with spaces, the text up to the first space is put in its own input_seq_name column for the BED (not to be confused w/ first name column, which specifies the pseudo-scaffold name). The remaining text in the input_seq_name after the first white space is split up at white space, delimited with given delimiter, and stored in a comment column of BED. Default = , (i.e. comma-delim).''')
parser.add_argument("-p", "--prefix", type=str, default=None, help='''Prefix for the output FASTA and BED files. Default: input FASTA file path/prefix up to the extension with .pseudo-scaffold.fasta/bed added.''')

args = parser.parse_args()

if args.prefix is None:
    args.prefix = '.'.join(args.fasta.split('.')[:-1]) + '.pseudo-scaffold'

# Open FASTA, build pseudo-scaffold, record seq coordinates, output new fasta and BED

scaf = ''
bed = ''
start = 0
with open(args.fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        desc = str(record.description).split()
        name = desc[0]
        if len(desc) > 1:
            comment = (args.commentdelim).join(desc[1:])
        else:
            comment = 'x' ## not using "-" to avoid confusion that it is a strand column
        seq = str(record.seq)
        seqlen = len(seq)
        end = start + seqlen
        scaf += seq
        bed += '\t'.join([str(e) for e in (args.name, start, end, name, comment, seqlen)]) + '\n'
        scaf += 'N'*args.length
        start = end + args.length

with open(args.prefix + '.fasta', 'w') as outfasta:
    outfasta.write('>'+args.name+'\n'+scaf.rstrip('N')+'\n')
with open(args.prefix + '.bed', 'w') as outbed:
    outbed.write(bed)

    
        
        
