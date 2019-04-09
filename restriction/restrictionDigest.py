#!/usr/bin/python
import sys
import argparse
import random
import os
from itertools import product

VersionInfo="""Version beta 2 
beta 1 version was created/updated on 11/14/2013 by John Urban.
beta 2 version was created/updated on 1/13/2019 by John Urban.
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
Optionally give the sequence a name by using the comma-separated format: name,sequence.
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
                    NOTE: the lengths output are calulated from restriction start site to start site.
                    However, the BED intervals will always extend the end from the start site to the restriction site end.
                    Therefore, if extracting sequences with these BED intervals, restriction site sequences will flank both ends.
                    For now --bed doesn't necessarily work for -t and -T -- even if you get an output.
                    Also, if you do get an output from -T, the BED intervals will only be start to start (i.e. ending restriction site not included).
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

parser.add_argument('--debug',
                    help="For development purposes, print debug/development messages to stderr.",
                    action="store_true", default=False)

parser.add_argument('--version', action='version', version=VersionInfo)

args = parser.parse_args()

##################################################################################################
'''---------------------------[ Debug functions ]----------------------------'''
##################################################################################################
def debug(msg):
    if args.debug:
        sys.stderr.write('[::restrictionDigest::] ' + msg +'\n')

##################################################################################################
'''---------------------------[ Early processing of arguments ]----------------------------'''
##################################################################################################
debug('Stage: Early processing of arguments.')
debug('Setting linear vs. circle.')
linear=True
if args.circular:
    linear = False

debug('Defining enzymes and restriction sites.')
enzymes = {'ApaI':'GGGCCC','BamHI':'GGATCC','EcoRI':'GAATTC','EcoRV':'GATATC','HindIII':'AAGCTT', 'KpnI':'GGTACC', 'PstI':'CTGCAG', 'XhoI':'CTCGAG', 'SalI':'GTCGAC', 'XbaI':'TCTAGA', 'BglI':'GCCNNNNNGGC', 'MstII':'CCTNAGG', 'HaeII':'RGCGCY'}
enzymes['Sau3AI'] = 'GATC'


if args.listEnzymes:
    print
    for key in sorted(enzymes.keys()):
        print key + "\t" + enzymes[key]
    print
    quit()



##########################################################################
'''---------------------------[ Functions ]----------------------------'''
##########################################################################
debug('Stage: Defining functions.')

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



##def allACGT(x):
##    xlen = len(x)
##    n = sum(e in 'ACGT' for e in x)
##    return n == xlen
##
def rseqify(rseq):
    ''' given seq return all possible seqs'''
    d = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'N':'ACGT'}
    return [''.join(j) for j in product(*(d[b] for b in rseq) ) ] 
    
def rseqify_list(x):
    ''' given list of restriction seqs, decode if needed for N,R,Y'''
    rseqs = []
    for rseq in x:
        rseqs += rseqify(rseq)
    return rseqs
    

    

## THIS WAS RETURNING 1-based starts....
##def substringPositions1(DNAstring, DNAsubstring):
##    """returns 1-based positions"""
##    window = len(DNAsubstring)
##    locations = []
##    for i in range(len(DNAstring)-window+1):
##        if DNAstring[i:i+window] == DNAsubstring:
##            locations = locations + [i+1] ## only diff from 0-based fxn
##    return locations


## NOW RETURNS 0-BASED STARTS
def substringPositions1(DNAstring, DNAsubstring):
    """returns 1-based positions"""
    window = len(DNAsubstring)
    locations = []
    for i in range(len(DNAstring)-window+1):
        if DNAstring[i:i+window] == DNAsubstring:
            locations = locations + [i] ## only diff from 0-based fxn
    return locations

def bedify_startslist(header, starts, motif):
    bedlist = []
    motiflen = len(motif)
    for start in starts:
        bedlist.append( (header, start, start+motiflen, motif) )
    return bedlist

def calculateLengths(starts, seq, lengths=[], beds=[], getbed=False, addXtoBedEnd=0):
    num = len(starts)
    if num == 1:
        lengths += [starts[0]]
        lengths += [len(seq)-starts[0]] 
        if getbed:
            beds.append( (0, starts[0]+addXtoBedEnd) )
            beds.append( (starts[0], len(seq) ) )
    if num > 1:
        lengths += [starts[0]]
        beds.append( (0, starts[0]) )
        for i in range(1, len(starts)):
            lengths += [starts[i] - starts[i-1]] ## Length is defined from start to start 
            if getbed:
                beds.append( (starts[i-1], starts[i]+addXtoBedEnd) )
        if len(starts) != starts[i-1]:
            lengths += [len(seq)-starts[i]] 
            if getbed:
                beds.append( (starts[i], len(seq)) )
    return lengths, beds



## Beta1 version
##def calculateCircularLengths(starts, seq, lengths=[], beds=[]):
##    """Need to be able to process plasmids"""
##    num = len(starts)
##    if num == 1:
##        lengths += [len(seq)]
##    if num > 1:
##        new1spot = starts[0]
##        seq = seq[new1spot :] + seq[:new1spot]
##        for i in range(len(starts)):
##            starts[i] = starts[i] - new1spot + 1
##        for i in range(1, len(starts)):
##            lengths += [starts[i] - starts[i-1]]
##        if len(starts) != starts[i-1]:
##            lengths += [len(seq)-starts[i]+1] ## I think +1 is wrong
##    return lengths


def calculateCircularLengths(starts, seq, lengths=[], beds=[], getbed=False, addXtoBedEnd=0):
    """Need to be able to process plasmids"""
    num = len(starts)
    if num == 1:
        lengths += [len(seq)]
        if getbed:
            beds.append( (starts[0], len(seq)+starts[0]+addXtoBedEnd) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
    if num > 1:
        for i in range(1, len(starts)): # skip first start: starts[0]
            lengths += [starts[i] - starts[i-1]]
            if getbed:
                beds.append( (starts[i-1], starts[i]+addXtoBedEnd) )
        # Last one = len-to-end + 0-to-start[0]
        lengths += [len(seq)-starts[i] + starts[0]] 
        if getbed:
            beds.append( (starts[i], len(seq)+starts[0]+addXtoBedEnd) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
    return lengths, beds

##def calculateCircularLengths(starts, seq, lengths=[], beds=[]):
##    """Need to be able to process plasmids"""
##    num = len(starts)
##    if num == 1:
##        lengths += [len(seq)]
##        beds.append( (starts[0], len(seq)+starts[0]) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
##    if num > 1:
##        new1spot = starts[0]
##        seq = seq[new1spot :] + seq[:new1spot]
##        for i in range(len(starts)):
##            starts[i] = starts[i] - new1spot + 1
##        for i in range(1, len(starts)):
##            lengths += [starts[i] - starts[i-1]]
##            beds.append( (starts[i-1], starts[i]) )
##        if len(starts) != starts[i-1]:
##            lengths += [len(seq)-starts[i]+1]
##    return lengths

def identity(x):
    return x

#"---------------------------[ Set up execution based on arguments ]----------------------------"

##########################################################################
'''---------------------------[ Execution ]----------------------------'''
##########################################################################
debug('Stage: Execution.')


##############################################################################################
'''---------------------------[ Loading sequence ]----------------------------'''
##############################################################################################
## Get Sequence
debug('Stage: Loading sequence.')
if args.fasta:
    debug('Setting sequence info from Fasta.')
    header, seq = fastaToPython(args.fasta)
    header = header[0].lstrip('>') ## for now, not doing multifastas
    seq = seq[0]
elif args.commandline:
    debug('Setting sequence info from command-line.')
    cl = args.commandline.split(',')
    if len(cl) == 1:
        header = 'command-line_sequence'
        seq = args.commandline
    elif len(cl) == 2:
        header = cl[0].lstrip('>')
        seq = cl[1]
    else:
        print "--commandline argument mis-used. Exiting..."
        quit()
seqlen = len(seq)
debug('Loaded ' + str(header) + ': ' + str(seqlen) + ' bp.')

##############################################################################################
'''---------------------------[ FIND RESTRICTION SITES ]----------------------------'''
##############################################################################################
debug('Stage: Finding restriction sites.')

## Process restriction site(s) as list
if args.restrictionEnzyme:
    args.restrictionSite = ','.join( rseqify_list([enzymes[e] for e in args.restrictionEnzyme.split(',')]) )
    

args.restrictionSite = args.restrictionSite.split(',') 



debug('Identifying restriction site starts in sequence: ' + str(header) )
if len(args.restrictionSite) > 1:
    debug('More than one enzyme detected - i.e. user provided a single enzyme or recognition sequence.')
    
    # startsList is a list of lists
    startsList = [substringPositions1(seq, args.restrictionSite[i]) for i in range(len(args.restrictionSite))]
    for e in startsList:
        debug('Detected ' + str(len(e)) + ' starts....')
    #bedlist = list of lists of tuples
    if args.bed:
        bedlist = [(header, startsList[i][j], startsList[i][j]+len(args.restrictionSite[i]), args.restrictionSite[i]) for j in range(len(startsList)) for i in range(len(args.restrictionSite))]


    if args.digestTogether: #Then combine list of lists into single list of starts and sort (as if only a single restrictionSite used)
        debug('DNA will be digested by more than one enzyme together in the same tube.')
        starts = sorted( [x for e in startsList for x in e] )
        debug('Detected ' + str(len(starts)) + ' sites overall.')

elif len(args.restrictionSite) == 1:
    debug('Only one enzyme detected - i.e. user provided more than one enzyme or recognition sequence.')

    # if only one restrictionSite, then just use a list of Starts
    starts = substringPositions1(seq, args.restrictionSite[0])
    debug('Detected ' + str(len(starts)) + ' sites overall.')




##############################################################################################
'''---------------------------[ CALCULATE LENGTHS ]----------------------------'''
##############################################################################################
debug('Stage: Calculate lengths.')


if args.digestTogether or len(args.restrictionSite) == 1: ## starts is a list
    debug('Single tube reaction.')
    if linear:
        debug('Calculating lengths given a linear molecule.')
        lengths, intervals = calculateLengths(starts, seq, lengths=[], getbed=args.bed, addXtoBedEnd=0)
    elif args.circular:
        debug('Calculating lengths given a circular molecule.')
        lengths, intervals = calculateCircularLengths(starts, seq, lengths=[], getbed=args.bed, addXtoBedEnd=0)

     
else: ## starts is a list of lists
    ## Get lengths list for each list in list of lists
    ## Store as list of lists
    debug('More than one reaction.')
    lengthsList = []
    intervals = []
    for i in range(len(startsList)):
        if linear:
            debug('Calculating lengths given a linear molecule for reaction ' + str(i) + '.')
            lengths, beds = calculateLengths(startsList[i], seq, lengths=[], getbed=args.bed, addXtoBedEnd=0)
            lengthsList += [lengths]
            intervals += [beds]
        elif args.circular:
            debug('Calculating lengths given a circular molecule for reaction ' + str(i) + '.')
            lengths, beds = calculateCircularLengths(startsList[i], seq, lengths=[], getbed=args.bed, addXtoBedEnd=0)
            lengthsList += [lengths]
            intervals += [beds]

    if args.loadTogether:
        debug('Loading restriction reactions into same well for gel electrophoresis.')
        ## then lengths is list of lengths you would see in a single lane
        ## combine the list of lists into a single list of lengths

        lengths = [x for e in lengthsList for x in e]
        #lengths = sorted( [x for e in startsList for x in e] )
    else:
        debug('NOT loading restriction reactions into same well for gel electrophoresis.')
        lengths = []


## ADD LENGTH OF UNDIGESTED MOLECULE?
if args.length and linear:
    debug('Intact linear molecule will be present.')
    lengths += [len(seq)]
elif args.length and args.circular:
    debug('Intact circular molecule(s) will be present.')
    lengths += ["super-coiled and nicked circles"]
    

## SORT LENGTHS INTO NEW LIST
debug('Sorting lengths.')
sortedLengths = sorted(lengths, reverse=True)







##############################################################################################
'''---------------------------[ OUTPUTS ]----------------------------'''
##############################################################################################
debug('Stage: Output.')

## For when >1 site given, but not to be put together
## Calculate separate reactions...
if args.bed:
    if len(args.restrictionSite) == 1:
        for interval in intervals:
            bed = [header, interval[0], interval[1], args.restrictionSite[0], interval[1]-interval[0], '+']
            print ('\t').join([str(e) for e in bed])
    else:
        debug('multiple-sites for BED still in dev. Need different solutions for digest together and sep... if doing sep, then just analyze one at a time... duh')
        quit()
##        for i in range(len(args.restrictionSite)):
##            for interval in intervals[i]:
##                bed = [header, interval[0], interval[1], args.restrictionSite[i]]
##                print ('\t').join([str(e) for e in bed])

elif len(args.restrictionSite) > 1 and not (args.loadTogether or args.digestTogether):
    debug('Lengths are in single well of gel.')
    for i in range(len(args.restrictionSite)):
        print args.restrictionSite[i]
        if args.proportion:
            debug('Calculating proportions for each fragment.')

            if args.length: ## This may seem redundant with section: ADD LENGTH OF UNDIGESTED MOLECULE? -- but its not b/c that is on "lengths" and this is using "lengthsList below
                print args.mass*len(seq)/float(len(seq))

            for j in range(len(lengthsList[i])):
                print args.mass*lengthsList[i][j]/float(len(seq))

            if i < len(args.restrictionSite)-1: 
                print ## Print space between current and next output.
        else:
            debug('Showing lengths of each fragment as seen in gel.')

            if args.length: ## This may seem redundant with section: ADD LENGTH OF UNDIGESTED MOLECULE? -- but its not b/c that is on "lengths" and this is using "lengthsList below
                print len(seq) 
                
            for j in range(len(lengthsList[i])):
                print lengthsList[i][j]

            if i < len(args.restrictionSite)-1: 
                print ## Print space between current and next output.
else:
    debug('Lengths are in separate wells of gel for different reactions.')
    if args.ordered:
        L = lengths
    else:
        L = sortedLengths
        
    for i in range(len(L)):
        if args.proportion:
            print args.mass*L[i]/float(len(seq))
        else:
            print L[i]


debug('All done.')

debug('DEVEL.')
quit()
if args.bed and len(args.restrictionSite) == 1: ## July 19, 2014 ##not done---
    debug('UNDEVELOPED BED CODE.')
    length = len(args.restrictionSite[0])
    center = len(args.restrictionSite[0])//2 ##center added to start
    ## if even, center is just chosen as 'right' from the left-right center pair
    if args.motifCenters:
        debug('UNDEVELOPED BED CODE: if motif centers.')
        for start in starts:
            bed = header[0][1:] + "\t" + str(start+center) + "\t" + str(start+center+1)
            sys.stdout.write(bed+"\n")
    else: ## report coordinates from start tp end of motif
        debug('UNDEVELOPED BED CODE: else (not motifcenters).')
        for start in starts:
            bed = header[0][1:] + "\t" + str(start) + "\t" + str(start+length)
            sys.stdout.write(bed+"\n")
    quit()

#'-----------------[ Clean up ]----------------------------------'
#if stdin:
#    os.remove(args.input)



#'----------------------[ Retired or developing ]-----------------------------------------'


