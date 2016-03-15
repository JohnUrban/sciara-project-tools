#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file and a file of entry names, return from the fastx file only those entries.

    """, formatter_class= argparse.RawTextHelpFormatter)


inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fasta', '-f',
                   type= str,
                   help='''Path to fasta file.''')
inputtype.add_argument('--stdin', action='store_true',
                   help='''Fastx is coming from stdin stream.''')

parser.add_argument('--quiet', '-q',
                   action='store_true',
                   help='''Suppress pct complete messages.''', default=False)

args = parser.parse_args()

############################################
if not args.fasta and not args.stdin:
    print "Specify --stdin or -f/--fasta file.fa"
    quit()


## functions
def find(string, char):
    return [i for i, ltr in enumerate(string) if ltr == char]

    
# Go through fasta file and return records that have IDs that match an element in set of names
if args.stdin:
    fastxFile = sys.stdin
else:
    fastxFile = open(args.fasta)
    
out = sys.stdout
msg = sys.stderr
nonredundant = {}
names = []

for record in SeqIO.parse(fastxFile, "fasta"):
    description=record.description
    seq=str(record.seq)
    if description in nonredundant.keys():
        if nonredundant[description] != seq:
            msg.write("Redundant name, but different seqs found: "+description+"\nIn these cases, sequences after first are ignored.\n\n")
    else:
        names.append(description)
        nonredundant[description] = seq

for name in names:
    out.write(">"+name+"\n"+nonredundant[name]+"\n")


fastxFile.close()
out.close()


