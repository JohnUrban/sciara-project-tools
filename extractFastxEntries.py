#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fasta/fastq file and a file of entry names, return from the fastx file only those entries (or only other entries).

Assumes each given name only occurs once in the file by default.

For very large fastx files, if you are looking to get >50% of the entries with supplied names,
it might be more efficient to supply the names of entries you do not want with "-e" option.

Similarly, for very large fastx files, if you are looking to exclude >50% of the entries with -e,
it might be more efficient to supply the names of entries you do want (without -e).

It more realistically depends on your workflow though.

"very large fastx files" are likely those with tens of millions of entries or more.
On a macbook pro, for a fasta file with 87509 scaffolds:
-- the approach of names to keep for 87256 names took 20 seconds
-- the approach of names to exclude for the complementary 253 entry names took 17 seconds
Not a big difference.

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

parser.add_argument('--exclude', '-e', action='store_true', default=False,
                   help='''All sequences in fastx file EXCEPT names given will be returned.''')

parser.add_argument('--multiple', '-m', action='store_true', default=False,
                    help=''' Only use this if given names possibly occur multiple times in the file. This is unusual.
By default, it is assumed that a given name only occurs once, so it is removed from the list when it is encountered.
This flag will turn that off.''')

parser.add_argument("--minlen",
                    type=int, default=0,
                    help='''Use this to extract reads >= int given. Default: 0.''')

parser.add_argument("--maxlen",
                    type=int, default=30000000000,
                    help='''Use this to extract reads <= int given. Default: 30 billion.''')

parser.add_argument('--separate', '-s',
                   action='store_true',
                   help='''Works with -n and -c. Puts each fasta record in its own fasta file called name.fasta (with given name). This overrides default behavior of printing to stdout.''', default=False)

parser.add_argument('-i', '--indexes', type=str, default=False,
                    help=''' Use this if you just want to extract specific entries according to their order of appearance.
For example, you may want only the first entry (use -i 0), or first 2 (-i 0:2), or the first 5 and last 5 (-i 0:5,N-5:N).
Or you can take every even numbered entry by -i 0:N:2, for example.
In the latter 2 cases, you need to know the values of N and N-5 for now.
In general, you can give comma-separated indexes and can give "slices" as colon-separated integers.
Indexes are 0-based like Python. Slices will go up to but not include end of range.''')

parser.add_argument('--head', type=int, default=False,
                    help=''' Use this if you just want to extract the head (first N bases) of all sequences.''')
##
##parser.add_argument('-U', '--touppercase', type=int, default=False,
##                    help=''' Will ensure all letters in sequence returned are uppercase.''')
##parser.add_argument('-L', '--tolowercase', type=int, default=False,
##                    help=''' Will ensure all letters in sequence returned are lowercase.''')

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
        if args.namesfile == "-" or args.namesfile == "stdin":
            nfile = sys.stdin
        else:
            nfile = open(args.namesfile)
        for line in nfile:
            names.add(line.rstrip())
    elif args.names:
        for name in args.names.split(","):
            names.add(name)
    setsize = len(names)
    pctdone = 0

    # Go through fasta file and return records that have IDs that match an element in set of names
    gatepct=10
    if not args.exclude:
        for record in SeqIO.parse(fastxFile, fastx):
            if record.id in names:
                if args.separate:
                    with open(record.id+".fasta", 'w') as f:
                        SeqIO.write(record, f, fastx)
                else:
                    SeqIO.write(record, out, fastx)
                if not args.multiple:
                    names.remove(record.id)
                pctdone += 100*1.0/setsize
                if pctdone >= gatepct and args.verbose:
                    msg.write(str(pctdone)+"% complete...."+jobname+"\n")
                    gatepct += 10
    else:
        for record in SeqIO.parse(fastxFile, fastx):
            if record.id not in names:
                if args.separate:
                    with open(record.id+".fasta", 'w') as f:
                        SeqIO.write(record, f, fastx)
                else:
                    SeqIO.write(record, out, fastx)
                pctdone += 100*1.0/setsize
                if pctdone >= gatepct and args.verbose:
                    msg.write(str(pctdone)+"% complete...."+jobname+"\n")
                    gatepct += 10
            else: #it is in names, so no need to check for this one anymore
                if not args.multiple:
                    names.remove(record.id)

elif args.indexes:
    extract = []
    indexes = args.indexes.split(",")
    for idx in indexes:
        idxrange = idx.split(":")
        if len(idxrange) == 2:
            extract += range(int(idxrange[0]), int(idxrange[1]))
        elif len(idxrange) == 3:
            extract += range(int(idxrange[0]), int(idxrange[1]), int(idxrange[2]))
        elif len(idxrange) == 1:
            extract.append(int(idxrange[0]))
    extract.sort()
    i = 0
    j = 0
    nidx = len(extract)
    for record in SeqIO.parse(fastxFile, fastx):
        if j < nidx and i == extract[j]:
            SeqIO.write(record, out, fastx)
            j += 1
        i += 1

elif args.head:
    for record in SeqIO.parse(fastxFile, fastx):
        record = record[:args.head]
        SeqIO.write(record, out, fastx)
    
elif args.minlen or args.maxlen:
    for record in SeqIO.parse(fastxFile, fastx):
        if len(record) >= args.minlen and len(record) <= args.maxlen:
            SeqIO.write(record, out, fastx)

fastxFile.close()
out.close()


