#!/usr/bin/env python2.7

import sys, argparse
from Bio import SeqIO
##from Bio.SeqIO import FastaIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION -
Input: a fastx file (usualy fasta).
Output: fastx file back w/ entries pasted together in order they appear.

Optional: Add NNNN-gaps of length X in between entries.


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
parser.add_argument('--gaplength', '-L',
                   type=int, default=0,
                   help='''Each fasta entry can be separated by NNN-gap of length X.
When given 0, entries are directly joined (no gaps). Default is 0 (no gapping).''')

parser.add_argument('--pastename', '-N',
                   type=str, default='pasted_seqs',
                   help='''The name of the pasted entry will include the names/descriptions of all entries in its comment section.
However, you can control what the main name of the pasted entries is -- called paste name.
By default the paste name is 'pasted_seqs'.
To use only the pastename (and exclude the collected info on all seq names), use --clean.
''')

parser.add_argument('--clean', 
                   action='store_true', default=False, 
                   help='''Use only paste name in output. See --pastename.''')

parser.add_argument('--outfastx', '-o',
                   type=str, required=True,
                   help='''Name of  file to store pasted fastx in. If you specify as 'stdout', it will print to stdout.
Cannot print both this and --outbed to stdout.''')

parser.add_argument('--outbed', '-O',
                   type=str, required=True,
                   help='''Name of  file to store BED locations of original sequences along pasted fastx.
If you specify as 'stdout', it will print to stdout.
Cannot print both this and out --outfastx to stdout.''')

args = parser.parse_args()

## ENSURE <= 1 STDOUT
try:
    if args.outfastx == 'stdout':
        assert args.outbed != 'stdout'
    elif args.outbed == 'stdout':
        assert args.outfastx != 'stdout'
except:
    print "Set either --outfasta or --outbed to stdout, but not both."
    quit()

## IF STDIN, REQUIRE FTYPE DECLARATION
if args.fastx == "-" or args.fastx == "stdin":
    args.fasta = sys.stdin
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
names = []
seqs = []
quals = []
starts = []
ends = []
lengths = []
currloc = 0 ## intitialize as currStart
for fa in SeqIO.parse(args.fastx , ftype):
    names.append( str(fa.description) )
    seqs.append( str(fa.seq) )
    starts.append(currloc) ## currloc is currStart
    seqlen = len(fa)
    currloc += seqlen ##adding seqlen to currStart gives currEnd
    ends.append(currloc)
    currloc += args.gaplength ## adding gapLen to currStart+seqlen sets currloc to nextstart
    if ftype == 'fastq':
        quals.append( str(fa.quals) )


    



## FORMAT FINAL NAME
if ftype == 'fasta':
    symbol = ">"
elif ftype == 'fastq':
    symbol = "@"
name = symbol + args.pastename
if not args.clean:
    name += '\t' + ('\t').join(names)

## FORMAT FINAL SEQ
seq = ('N'*args.gaplength).join(seqs)

## IF FQ, FORMAT FINAL QUAL STRING
if ftype == 'fastq':
    ## ! should encode Q=0
    qual = ('!'*args.gaplength).join(quals)


## WRITING FILES
## DETERMINE OUTPUT LOCATIONS FOR OUTFASTX 
if args.outfastx == 'stdout':
    fout = sys.stdout
else:
    fout = open(args.outfastx, 'w')
## OUTPUT FASTX
fout.write( name + "\n") ##1st line
fout.write( seq + "\n") ## 2nd line
if ftype == 'fastq':
    fout.write( '+\n' ) ##3rd line
    fout.write(qual + "\n") ##4th line
fout.close()


## DETERMINE OUTPUT LOCATIONS FOR OUTBED 
if args.outbed == 'stdout':
    bout = sys.stdout
else:
    bout = open(args.outbed, 'w')
## OUTPUT BED
for i in range(len(starts)):
    start = starts[i]
    end = ends[i]
    desc = ('_').join(names[i].split())
    entry = ('\t').join([str(e) for e in [args.pastename, start, end, desc]])
    bout.write(entry + "\n")
bout.close()


