#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file and a file of entry names, return from the fastx file only those entries.

    """, formatter_class= argparse.RawTextHelpFormatter)


inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fastx', '-f',
                   type= str,
                   help='''Path to fasta/fastq file.''')
inputtype.add_argument('--stdin', action='store_true',
                   help='''Fastx is coming from stdin stream.''')


filetype = parser.add_mutually_exclusive_group()
filetype.add_argument('--fa', action='store_true',
                   help='''Explicitly define as fasta input.''')
filetype.add_argument('--fq', action='store_true',
                   help='''Explicitly define as fastq input.''')


parser.add_argument('--verbose', '-v',
                   action='store_true',
                   help='''Provides jobname and pct complete messages while running (only applicable to extracting entries based on names).''', default=False)

names = parser.add_mutually_exclusive_group()

names.add_argument('--namesfile', '-n', type=str,
                   help='''Path to file with fasta/fastq entry names -- one name per line''')

names.add_argument('--names', '-c', type=str, help='Enter comma-separated names at command line with this flag')


parser.add_argument("--minlen",
                    type=int, default=0,
                    help='''Use this to extract reads >= int given. Default: 0.''')

parser.add_argument("--maxlen",
                    type=int, default=30000000000,
                    help='''Use this to extract reads <= int given. Default: 30 billion.''')

parser.add_argument('--separate', '-s',
                   action='store_true',
                   help='''Works with -n and -c. Puts each fasta record in its own fasta file called name.fasta (with given name). This overrides default behavior of printing to stdout.''', default=False)


args = parser.parse_args()

############################################
'''           functions                '''
############################################

def find(string, char):
    return [i for i, ltr in enumerate(string) if ltr == char]

############################################
'''           check options              '''
############################################

## ensure that one of these are present
## quit if no -n or -c
    ##if not args.namesfile and not args.names:
    ##    print "-n or -c is required. Use -h for help."
    ##    quit()

## ensure one of these is used, and if stdin then also explicity specify filetype    
if not args.fastx and not args.stdin:
    print "Specify --stdin or -f/--fastx file.fx"
    quit()
if args.stdin and not (args.fa or args.fq):
    print "When using --stdin, filetype should be explicitly defined --fa or --fq."
    quit()

# file or stdin?
if args.stdin:
    fastxFile = sys.stdin
else:
    fastxFile = open(args.fastx)

# FASTA or FASTQ?
if args.fa:
    fastx = "fasta"
elif args.fq:
    fastx = "fastq"
elif args.fastx: ## can figure out from file (right now not from stdin due to non-seekability)
    line1 = fastxFile.next()[0]
    if line1[0] == ">":
        fastx = "fasta"
    elif line1[0] == "@":
        fastx = "fastq"
    fastxFile.seek(0)
else:
    print "Expected fasta or fastq. File given needs to be reformatted if user thinks it is."
    quit()



# job name
if args.verbose:
    pathInd = []
    if args.fastx:
        pathInd = find(args.fastx, "/")
    if len(pathInd) == 0:
        if args.stdin:
            jobname = "stdin"
        else:
            jobname = args.fastx
    else:
        jobname = args.fastx[max(pathInd)+1:]
    jobname += "+"
    if args.namesfile:
        pathInd = find(args.namesfile, "/")
        if len(pathInd) == 0:
            jobname += args.namesfile
        else:
            jobname += args.namesfile[max(pathInd)+1:]
    elif args.names:
        jobname += "commandline"




############################################
'''           execute                '''
############################################    

out = sys.stdout
msg = sys.stderr

if args.namesfile or args.names:
    ## Make set of record IDs (names)
    names = set()
    if args.namesfile:
        for line in open(args.namesfile):
            names.add(line.rstrip())
    elif args.names:
        for name in args.names.split(","):
            names.add(name)
    setsize = len(names)
    pctdone = 0

    # Go through fasta file and return records that have IDs that match an element in set of names
    gatepct=10
    for record in SeqIO.parse(fastxFile, fastx):
        if record.id in names:
            if args.separate:
                with open(record.id+".fasta", 'w') as f:
                    SeqIO.write(record, f, fastx)
            else:
                SeqIO.write(record, out, fastx)
            names.remove(record.id)
            pctdone += 100*1.0/setsize
            if pctdone >= gatepct and args.verbose:
                msg.write(str(pctdone)+"% complete...."+jobname+"\n")
                gatepct += 10

elif args.minlen or args.maxlen:
    for record in SeqIO.parse(fastxFile, fastx):
        if len(record) >= args.minlen and len(record) <= args.maxlen:
            SeqIO.write(record, out, fastx)

fastxFile.close()
out.close()


