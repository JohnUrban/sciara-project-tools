#!/usr/bin/python
import sys
import argparse
import random
import os
VersionInfo="""Version beta 1 
This version was created/updated on 11/14/2013 by John Urban.
 """

parser = argparse.ArgumentParser(description=VersionInfo+"""

DESCRIPTION

    Provide a fasta file or sequence at command line and a restriction site sequence (or other sub-sequence).
    This will output approximations of expected band sizes by cutting approximately in the middle of sub-sequence
    
EXAMPLE INPUT:
        
        
EXAMPLE OUTPUT:
    


Example usages:
Find lengths of HindIII sites if digest Lambda genome

restrictionDigest.py -f lambda.fa -r AAGCTT

23130
9416
6682
4362
2322
2027
564

restrictionDigest.py -f lambda.fa -r AAGCTT -l
48502
23130
9416
6682
4362
2322
2027
564




NOTES:
    This should only be used on small sequences.
    Bowtie or Bowtie2 should be used to map restriction sites genome-wide.

    See the restriction sites script I ALREADY MADE before continuing this.


TODO:


        

    """,
    formatter_class= argparse.RawTextHelpFormatter)

inputSeq = parser.add_mutually_exclusive_group()
inputSeq.add_argument('-c', '--commandline',
                    type=str,
                    help='''Entering sequence at command line.
Can also pipe in the output -- use "-" or "stdin" instead.
                    ''')

inputSeq.add_argument('-f', '--fasta',
                    type=str,
                    help='''Providing a sequence in fasta format.
Can also pipe in the output -- use "-" or "stdin" instead.
                    ''')
restriction = parser.add_mutually_exclusive_group()
restriction.add_argument('-r', '--restrictionSite',
                    type=str,
                    help='''The restriction site(s) or sub-sequence(s).
                    Either put 1 subsequence, or a comma-separated list.
                    If giving multiple, see -t to print out together instead of individually.
                    See -T to digest together.''')

restriction.add_argument('-R', '--restrictionEnzyme',
                    type=str,
                    help='''The restriction enzyme(s).
                    Either put 1 enzyme, or a comma-separated list.
                    If giving multiple, see -t to print out together instead of individually.
                    See -T to digest together.''')

restriction.add_argument('-L', '--listEnzymes',
                         help='List available enzymes.',
                         action='store_true')

parser.add_argument('-t', '--loadTogether',
                    help='''If multiple restriction sites given, this flag says to digest separately, then combine in one tube.
                    What would that ladder/fingerprint look like.''',
                    action="store_true")

parser.add_argument('-T', '--digestTogether',
                    help='''If multiple restriction sites given, this flag says to digest separately, then combine in one tube.
                    What would that ladder/fingerprint look like.''',
                    action="store_true")

parser.add_argument('-l', '--length',
                    help='''Also include length of original sequence in the output.''',
                    action="store_true")

parser.add_argument('-P', '--proportion',
                    help='''Report band lengths as proportions of full length.''',
                    action="store_true")

parser.add_argument('-M', '--mass',
                    type=float,
                    help='''Ng, Ug (or any unit) amount to multiply proportion against.
                    Only works when -P is specified also.
                    Default is 1.
                    Can use 100 to convert proportions to percents.
                    Can use 500 (e.g. if digesting 500 ng) to see how many ng each band will have.''',
                    default=1.0)

parser.add_argument('-C', '--circular',
                    help='''USe this flag if the input sequence should be considered as a circle, not linear.''',
                    action="store_true")

parser.add_argument('-B', '--bed',
                    help='''Tells it to output as a BED file.
                    ''',
                    action="store_true")

parser.add_argument('-m', '--motifCenters',
                    help='''Tells it to output as a BED file.
                    ''',
                    action="store_true")

parser.add_argument('-S', '--sam',
                    help='''Tells it that output as a SAM file.
                    ''',
                    action="store_true")

parser.add_argument('-o', '--outName',
                    type=str,
                    help='''The name of the output file (optional).
Default is stdout.
                    ''',
                    default="stdout")

parser.add_argument('-O', '--ordered',
                    action='store_true',
                    default=False,
                    help='''Return lengths in order encountered.
                    ''')

parser.add_argument('-v', '--verbose',
                    help="Use to increase verbosity while program is running",
                    action="store_true")
parser.add_argument('--version', action='version', version=VersionInfo)

args = parser.parse_args()


#"---------------------------[ Check and pare arguments ]----------------------------"
linear=True
if args.circular:
    linear = False

if args.listEnzymes:
    print
    print 'ApaI         GGGCCC'
    print 'BamHI        GGATCC'
    print 'EcoRI        GAATTC'
    print 'EcoRV        GATATC'
    print 'HindIII      AAGCTT'
    print 'KpnI         GGTACC'
    print
    quit()
enzymes = {'ApaI':'GGGCCC','BamHI':'GGATCC','EcoRI':'GAATTC','EcoRV':'GATATC','HindIII':'AAGCTT', 'KpnI':'GGTACC'}

if args.restrictionEnzyme:
    args.restrictionSite = enzymes[args.restrictionEnzyme]
#"---------------------------[ Functions ]----------------------------"

def fastaToPython(fastaFile):
    seqs = open(fastaFile, 'r')
    fasta = seqs.readlines()
    seqs.close()
    headers = []
    DNAstrings = []
    index=0
    ## puts names in 'headers'
    ## puts sequence in 'DNAstrings'
    ##      with same index as corresponding name in 'headers'
    ## To do this, I combine all lines in between '>' header lines
    ## Moreover, I make into a single string not separated by '\n'
    for element in fasta:
        if element[0] == '>':
            if element[-1] == '\n':
                headers = headers + [element[:-1]]
            else:
                headers = headers + [element]
            DNAstrings = DNAstrings + ['']
            index+=1
        else:
            if element[-1] == '\n':
                DNAstrings[index-1] = DNAstrings[index-1] + element[:-1]
            else:
                DNAstrings[index-1] = DNAstrings[index-1] + element
    
    return headers, DNAstrings

def substringPositions1(DNAstring, DNAsubstring):
    """returns 1-based positions"""
    window = len(DNAsubstring)
    locations = []
    for i in range(len(DNAstring)-window+1):
        if DNAstring[i:i+window] == DNAsubstring:
            locations = locations + [i+1] ## only diff from 0-based fxn
    return locations


def calculateLengths(starts, seq, lengths=[]):
    num = len(starts)
    if num == 1:
        lengths += [starts[0]]
        lengths += [len(seq)-starts[0]]
    if num > 1:
        lengths += [starts[0]]
        for i in range(1, len(starts)):
            lengths += [starts[i] - starts[i-1]]
        if len(starts) != starts[i-1]:
            lengths += [len(seq)-starts[i]+1]
    return lengths

def calculateCircularLengths(starts, seq, lengths=[]):
    """Need to be able to process plasmids"""
    num = len(starts)
    if num == 1:
        lengths += [len(seq)]
    if num > 1:
        new1spot = starts[0]
        seq = seq[new1spot :] + seq[:new1spot]
        for i in range(len(starts)):
            starts[i] = starts[i] - new1spot + 1
        for i in range(1, len(starts)):
            lengths += [starts[i] - starts[i-1]]
        if len(starts) != starts[i-1]:
            lengths += [len(seq)-starts[i]+1]
    return lengths

#"---------------------------[ Set up execution based on arguments ]----------------------------"

#'-----------------[ execution ]----------------------------------'

##Process restriction site(s)
args.restrictionSite = args.restrictionSite.split(',') ## list

## Get Sequence
if args.fasta:
    header, seq = fastaToPython(args.fasta)
    seq = seq[0]
elif args.commandline:
    seq = args.commandline


if len(args.restrictionSite) > 1:
    startsList = []
    for i in range(len(args.restrictionSite)):
        startsList += [substringPositions1(seq, args.restrictionSite[i])]
    if args.digestTogether:
        starts = []
        for i in range(len(startsList)):
            starts += startsList[i]
        starts.sort()
elif len(args.restrictionSite) == 1:
    starts = substringPositions1(seq, args.restrictionSite[0])

if args.bed and len(args.restrictionSite) == 1: ## July 19, 2014 ##not done---
    length = len(args.restrictionSite[0])
    center = len(args.restrictionSite[0])//2 ##center added to start
    ## if even, center is just chosen as 'right' from the left-right center pair
    if args.motifCenters:
        for start in starts:
            bed = header[0][1:] + "\t" + str(start+center) + "\t" + str(start+center+1)
            sys.stdout.write(bed+"\n")
    else: ## report coordinates from start tp end of motif
        for start in starts:
            bed = header[0][1:] + "\t" + str(start) + "\t" + str(start+length)
            sys.stdout.write(bed+"\n")
    quit()

### FIND LENGTHS
if args.digestTogether or len(args.restrictionSite) == 1:
    if linear:
        lengths = calculateLengths(starts, seq, lengths=[])
    elif args.circular:
        lengths = calculateCircularLengths(starts, seq, lengths=[])
        
else:
    lengthsList = []
    for i in range(len(startsList)):
        if linear:
            lengthsList += [calculateLengths(startsList[i], seq, lengths=[])]
        elif args.circular:
            lengthsList += [calculateCircularLengths(startsList[i], seq, lengths=[])]
    if args.loadTogether:
        lengths = []
        for i in range(len(lengthsList)):
            lengths += lengthsList[i]
    else:
        lengths = []


if args.length and linear:
    lengths += [len(seq)]
elif args.length and args.circular:
    lengths += ["super-coiled and nicked circles"]
    

sortedLengths = sorted(lengths, reverse=True)

## For when >1 site given, but not to be put together
if len(args.restrictionSite) > 1 and not (args.loadTogether or args.digestTogether):
    for i in range(len(args.restrictionSite)):
        print args.restrictionSite[i]
        if args.proportion:
            if args.length:
                print args.mass*len(seq)/float(len(seq))
            for j in range(len(lengthsList[i])):
                print args.mass*lengthsList[i][j]/float(len(seq))
            if i < len(args.restrictionSite)-1: ## TOFINISH??
                print             
        else:
            if args.length:
                print len(seq)
            for j in range(len(lengthsList[i])):
                print lengthsList[i][j]
            if i < len(args.restrictionSite)-1: ## TOFINISH??
                print 
else:
    if args.ordered:
        L = lengths
    else:
        L = sortedLengths
        
    for i in range(len(L)):
        if args.proportion:
            print args.mass*L[i]/float(len(seq))
        else:
            print L[i]


#'-----------------[ Clean up ]----------------------------------'
#if stdin:
#    os.remove(args.input)



#'----------------------[ Retired or developing ]-----------------------------------------'


