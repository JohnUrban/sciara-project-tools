#!/usr/bin/env python2.7

import sys, argparse
from Bio import SeqIO
##from Bio.SeqIO import FastaIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION -
Input:
    1. a fastx file (usualy fasta).
    2. file containing desired order of fastx entries (1 name per line)
Output: fastx file back w/ entries ordered according to desired order.

This reads in all fastx entries, so be careful w/ large files w/ many entries.

For now this requires the entire name/description for each entry as it appears in the fastx file.
So that is everything MINUS the ">" or "@" symbol.
Example:
>this is a name
ACGT

...should appear in order file as "this is a name" as opposed to "this" (for example).

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fastx', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file in FASTX format.
Declaring input as fasta or fastq takes precedence (--fasta, --fastq).
If not declared, Fasta vs. Fastq can be auto-detected by file extension.
Fasta: .fasta, .fa, .fna, .fpa
Fastq: .fastq, .fq
If that does not work, it attempts auto-detection by looking at the first character of the first line.
">" = fasta
"@" = fastq
''')
ftype = parser.add_mutually_exclusive_group()
ftype.add_argument('--fasta', '-fa',
                   action='store_true', default=False, 
                   help='''Declare input as FASTA format.''')
ftype.add_argument('--fastq', '-fq',
                   action='store_true', default=False, 
                   help='''Declare input as FASTQ format.''')

parser.add_argument('--outfastx', '-o',
                   type=str, default='stdout',
                   help='''Name of  file to store pasted fastx in. If you specify as 'stdout', it will print to stdout.
Defaults to stdout.''')

parser.add_argument('--orderfile', '-O',
                   type=str, required=True,
                   help='''Name of  file to that has desired order of fastx entries.''')


parser.add_argument('--qualoffset', '-Q', type=int, default=33,
                    help='''Fastq qual strings are read in assuming Phred+33. This option only affects output for now.
                    So it should be kept as default of 33.''')


args = parser.parse_args()


## IF STDIN, REQUIRE FTYPE DECLARATION
if args.fastx == "-" or args.fastx == "stdin":
    args.fastx = sys.stdin
    try:
        assert args.fasta or args.fastq
    except:
        print "For stdin as input, you need declare the file type with --fasta or --fastq"
        quit()

## DETERMINE FTYPE
ext = args.fastx.strip().split('.')[-1]
if args.fasta or (ext in ('fasta', 'fa', 'fna', 'fpa')):
    ftype = 'fasta'
elif args.fastq or (ext in ('fastq', 'fq')):
    ftype = 'fastq'
else:
    ftype_clue = open(args.fastx, 'r')
    symbol = ftype_clue.readline()[0]
    if symbol == ">":
        ftype = 'fasta'
    elif symbol == "@":
        ftype = 'fastq'
    else:
        print "Could not determine file type. Try declaring file type w/ --fasta or --fastq."
        quit()


    

## COLLECT SEQ INFO AND DEFINE LOCATIONS OF INDIV SEQS
seqs = {}
for fa in SeqIO.parse(args.fastx , ftype):
    seqs[str(fa.description)] = [str(fa.seq)]
    if ftype == 'fastq':
        qualstr = ('').join([str(unichr(q + args.qualoffset)) for q in fa.letter_annotations['phred_quality']])
        seqs[str(fa.description)] += ['+', qualstr]


    

## WRITING FILES
## DETERMINE OUTPUT LOCATIONS FOR OUTFASTX 
if args.outfastx == 'stdout':
    fout = sys.stdout
else:
    fout = open(args.outfastx, 'w')

## OUT SYMBOL
if ftype == 'fasta':
    symbol = ">"
elif ftype == 'fastq':
    symbol = "@"

## OUTPUT FASTX IN DESIRED ORDER
with open(args.orderfile, 'r') as ordfile:
    for line in ordfile:
        name = line.strip()
        fout.write( symbol + name + "\n" )
        fout.write( ('\n').join(seqs[name]) + '\n' )


