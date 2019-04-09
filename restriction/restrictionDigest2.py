#!/usr/bin/env python
import sys, re, string
import argparse
import random
import os
from itertools import product
from Bio import SeqIO
from collections import defaultdict

VersionInfo="""Version beta 2.1
beta 1 version was created/updated on 11/14/2013 by John Urban.
beta 2 version was created/updated on 1/13/2019 by John Urban.
beta 2.1 version was created/updated on 1/17/2019 by John Urban
    An overhaul of code is planned, so it was renamed to restrictionDigest2.py.
    Plans include:
        - Make use of regex module
            - only weakness would be non-overlapping starts
            - however, that is not typical for r sites and it would ignore tiny fragments not seen on gel
        - Treat multiple enzymes and single enzymes in same way
            - i.e. if use a list of lists for > 1, use it also for 1
        - Use dictionaries to:
            - keep enzymes associated with restriction seqs
            - restriction seqs associated with enzymes being used
                - note that subsetting is important since enzymes can have same seqs
            - keep start site info associated to the enzyme or seq
        - Either use SeqIO from Bio or my own fasta class to parse FASTA file
        - Be able to do this on multiple sequences


Version 2 End -- I now am on Version 3.
Summary:
I did end up making a lot of changes.
Any pipeline that originally ran restrictionDigest.py should either still use that or figure out the differences needed for Ver2.

I did not end up implementing the things listed in the overhaul exactly as planned.
Nonetheless, here are updates in Ver2:
- outputs BED files
- outputs rmap files
- works on multi-fastas
- default output now includes a header line (seqname)
- more verbose/debug messages
- does not load sequences all into memory any more; uses SeqIO instead and processes one at a time
- automatically updates enzymes dictionary with "rseqified"/expanded sequences (for those w/ N,R,Y in seq).
- automatically checks if revcomp is same or not (adds it if not)
- initial regex functions added

Version 3 plans:
I will begin to expand upon the regex approach.
- The regex pipeline will be optional
- Regexes will be used with re module for starts
- Ultimately, all recognition sequences will be converted into a single regex to search sequence once.
    - Be careful if recog. seqs have overlapping prefixes/suffixes.
- 

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
                    
                    ''',
                    action="store_true")


parser.add_argument('-rmap', '--rmap',
                    help='''Tells it to output as a rmap file: SeqName Length NumFrags Fraglen1 Fraglen2 .... FraglenN.
                    This outputs length-ordered maps (as seen on a gel) by default as per the rest of the program.
                    To see as a sequence-ordered map as for optical maps, use -O/--ordered.
                    
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
                    help='''Use to increase verbosity while program is running.
                    Not used at the moment. Use --debug instead.''',
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
    if args.debug or args.verbose:
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
enzymes['HpaII'] = 'CCGG'
enzymes['MspI'] = 'CCGG'

## I need this here until I implement regular expressions
enzymes['BglI_1'] = 'GCCAATTCGGC'
enzymes['BglI_2'] = 'GCCCAGTGGGC' #GCCC AGC TGGC
enzymes['BglI_3'] = 'GCCCAGCTGGC'

enzymes['BssSI'] = 'CACGAG'
enzymes['BssSIrc'] = 'CTCGTG'

enzymes['BspQI'] = 'GCTCTTC'
enzymes['BspQIrc'] = 'GAAGAGC'


## These below are needed for output stage (BED)
enzymes['5p'] = '*'
enzymes['3p'] = '*'

enzymelookup = defaultdict(list) ## THIS IS FILLED OUT ON THE FLY - I've added this in anticipation of the Regex approach of searching a sequence in one pass




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


##def fastaToPython(fastaFile):
##    seqs = open(fastaFile, 'r')
##    fasta = seqs.readlines()
##    seqs.close()
##    headers = []
##    DNAstrings = []
##    index=0
##    ## puts names in 'headers'
##    ## puts sequence in 'DNAstrings'
##    ##      with same index as corresponding name in 'headers'
##    ## To do this, I combine all lines in between '>' header lines
##    ## Moreover, I make into a single string not separated by '\n'
##    for element in fasta:
##        if element[0] == '>':
##            if element[-1] == '\n':
##                headers = headers + [element[:-1]]
##            else:
##                headers = headers + [element]
##            DNAstrings = DNAstrings + ['']
##            index+=1
##        else:
##            if element[-1] == '\n':
##                DNAstrings[index-1] = DNAstrings[index-1] + element[:-1]
##            else:
##                DNAstrings[index-1] = DNAstrings[index-1] + element
##    
##    return headers, DNAstrings


##def fastaToPython(fastaFile, seqdict={}):
##    with open(fastaFile, 'r') as seqs:
##        fasta = seqs.readlines()
##    ## puts names in 'headers'
##    ## puts sequence in 'DNAstrings'
##    ##      with same index as corresponding name in 'headers'
##    ## To do this, I combine all lines in between '>' header lines
##    ## Moreover, I make into a single string not separated by '\n'
##    for element in fasta:
##        if element[0] == '>':
##            header = element.strip()
##            debug("Loading : " + header )
##            seqdict[header] = ''
##        else:
##            seqdict[header] += element.strip()
##    return seqdict

def fastaToPython(fastaFile, seqdict={}):
    with open(fastaFile, 'r') as fasta:
        for element in fasta:
            if element[0] == '>':
                header = element.strip()
                debug("Loading : " + header )
                seqdict[header] = ''
            else:
                seqdict[header] += element.strip()
    return seqdict

##def allACGT(x):
##    xlen = len(x)
##    n = sum(e in 'ACGT' for e in x)
##    return n == xlen
##

def rseqify(rseq):
    ''' given seq return all possible seqs'''
    d = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'N':'ACGT'}
    return [''.join(j) for j in product(*(d[b] for b in rseq) ) ] 


def regexify(rseq):
    ''' given seq return all possible seqs'''
    d = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'N':'ACGT'}
    regex = ''
    oldb = rseq[0]
    cnt = 1
    for i in range(1,len(rseq)):
        b = rseq[i]
        if b == oldb:
            cnt += 1
            continue

        
        if len(d[oldb]) == 1:
            regex += oldb
        else:
            regex += '[' + d[oldb] + ']'

        if cnt > 1:
            regex += '{' + str(cnt) + '}'
        cnt = 1
        oldb = b

    ## PROCESS LAST ONE
    if len(d[b]) == 1:
        regex += b
    else:
        regex += '[' + d[b] + ']'

    if cnt > 1:
        regex += '{' + str(cnt) + '}'

    return regex


def revcomp(rseq):
    intab=  'actguACTGUNRY'
    outtab= 'tgacaTGACANYR'
    transtab = string.maketrans(intab, outtab)
    return rseq.translate(transtab)[-1::-1]    

##def regexifyrevcomp(rseq):
##    rc = revcomp(rseq)
##    if rc == rseq:
##        return regexify(rseq)
##    return regexify(rc)

def rseqify_list(x):
    ''' given list of restriction seqs, decode if needed for N,R,Y'''
    rseqs = []
    for rseq in x:
        rseqs += rseqify(rseq)
    return rseqs
    

def rseqify_list2(x):
    ''' given list of restriction enzymes, decode if needed for N,R,Y'''
    enames = []
    rseqs = []
    for e in x:
        seqs = rseqify(enzymes[e])
        lenseqs = len(seqs)
        if lenseqs == 1:
            enames.append( e )
        elif lenseqs > 1:
            enames += [e+"_"+str(i) for i in range(lenseqs)]
        rseqs += seqs

        #RCs
        rc = revcomp(enzymes[e])
        if rc != enzymes[e]:
            rcseqs = rseqify(rc)
            lenrcseqs = len(rcseqs)
            if lenrcseqs == 1:
                enames.append( e + '_rc' )
            elif lenrcseqs > 1:
                enames += [e+"_rc_"+str(i) for i in range(lenrcseqs)]
            rseqs += rcseqs
                
            
    return enames, rseqs    


def regexify_list(x):
    ''' given list of restriction enzymes, decode if needed for N,R,Y'''
    enames = []
    rseqs = []
    for e in x:
        enames.append( e )
        rseqs,append( regexify(enzymes[e]) )

        #RCs
        rc = revcomp(enzymes[e])
        if rc != enzymes[e]:
            enames.append( e + '_rc' )
            rseqs.append( regexify(rc) )        
    return enames, rseqs    


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
    """returns 0-based positions"""
    window = len(DNAsubstring)
    locations = []
    for i in range(len(DNAstring)-window+1):
        if DNAstring[i:i+window] == DNAsubstring:
            locations = locations + [i] ## only diff from 0-based fxn
    return locations


def regexStartPositions(DNAstring, DNAsubstring):
    ''' Has same outout as 0-based substringPositions1 fxn'''
    regex = re.compile(DNAsubstring)
    locations = []
    for m in re.finditer(regex, DNAstring):
        locations.append(m.start())
    return locations





def seqlist2regexlist(l):
    ''' given list of sequences, return list of regexes'''
    return 

def regexlist2regex(l):
    ''' given list of sequences, return list of regexes'''
    return '|'.join(l)

def bedify_startslist(header, starts, motif):
    bedlist = []
    motiflen = len(motif)
    for start in starts:
        bedlist.append( (header, start, start+motiflen, motif) )
    return bedlist

def calculateLengths(starts, seq, lengths=[], beds=[], getbed=False, addXtoBedEnd=[], labels=[]):
    num = len(starts)
    #debug("calculateLengths Fxn: "+str(num)+" "+str(len(lengths))+" "+str(len(beds)))
    if num == 1:
        #debug("calculateLengths Fxn: 1")
        lengths += [starts[0]]
        lengths += [len(seq)-starts[0]] 
        if getbed:
            beds.append( (0, starts[0]+addXtoBedEnd[0], '5p:'+labels[0]) )
            beds.append( (starts[0], len(seq), labels[0]+':3p' ) )
    if num > 1:
        #debug("calculateLengths Fxn: MORE THAN 1")
        lengths += [starts[0]]
        beds.append( (0, starts[0]+addXtoBedEnd[0], '5p:'+labels[0]) )
        for i in range(1, len(starts)):
            lengths += [starts[i] - starts[i-1]] ## Length is defined from start to start 
            if getbed:
                beds.append( (starts[i-1], starts[i]+addXtoBedEnd[i], labels[i-1]+':'+labels[i]) )
        if len(starts) != starts[i-1]:
            lengths += [len(seq)-starts[i]] 
            if getbed:
                beds.append( (starts[i], len(seq), labels[i]+':3p') )
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


def calculateCircularLengths(starts, seq, lengths=[], beds=[], getbed=False, addXtoBedEnd=[], labels=[]):
    """Need to be able to process plasmids"""
    num = len(starts)
    if num == 1:
        lengths += [len(seq)]
        if getbed:
            beds.append( (starts[0], len(seq)+starts[0]+addXtoBedEnd[0], labels[0]+':'+labels[0]) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
    if num > 1:
        for i in range(1, len(starts)): # skip first start: starts[0]
            lengths += [starts[i] - starts[i-1]]
            if getbed:
                beds.append( (starts[i-1], starts[i]+addXtoBedEnd[i], labels[i-1]+':'+labels[i]) )
        # Last one = len-to-end + 0-to-start[0]
        lengths += [len(seq)-starts[i] + starts[0]] 
        if getbed:
            beds.append( (starts[i], len(seq)+starts[0]+addXtoBedEnd[i] , labels[i]+':'+labels[0]) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
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


class Digest(object):
    def __init__(self, header, seq, enzymenames, recseqs):
        '''
        header is single fasta header or other string
        seq is single associated fasta sequence string
        names is a list of N strings for enzyme names or naming of recseqs somehow
        recseqs is a list of N associated strings/recognition sequences (can do single enzymedict in future, but code is made for comma-sep list atm)
        
        '''
        self.header = seqname
        self.seq = seq
        self.names = enzymenames
        self.recseqs = enzymeseqs
        self.n = len(enzymenames)
        assert self.n == len(enzymeseqs)
        self.niter = range(self.n)
        # Get Start positions of each enzymeseq
        self.startsList = [substringPositions1(self.seq, self.recseqs[i]) for i in self.niter]
        





#"---------------------------[ Set up execution based on arguments ]----------------------------"

##########################################################################
'''---------------------------[ Execution ]----------------------------'''
##########################################################################
debug('Stage: Execution.')

################################################################################################
##'''---------------------------[ Loading sequence ]----------------------------'''
################################################################################################
#### Get Sequence
##debug('Stage: Loading sequence(s).')
##if args.fasta:
##    debug('Processing sequence info from Fasta.')
##    #headers, seqs = fastaToPython(args.fasta)
##    headers, seqs = map(list, zip(*fastaToPython(args.fasta).iteritems()))
####    header = header[0].lstrip('>') ## for now, not doing multifastas
####    seq = seq[0]
##elif args.commandline:
##    debug('Processing sequence info from command-line.')
##    cl = args.commandline.split(',')
##    if len(cl) == 1:
##        headers = ['command-line_sequence']
##        seqs = [args.commandline]
##    elif len(cl) == 2:
##        header = cl[0].lstrip('>')
##        seq = cl[1]
##    else:
##        print "--commandline argument mis-used. Exiting..."
##        quit()
##
##nseq=len(seqs)
##nhead = len(headers)
##debug('Loaded ' + str(nhead) + ' headers and ' + str(nseq) + ' associated sequences.')
##assert nseq == nhead



##############################################################################################
'''---------------------------[ PROCESS RESTRICTION SITES ]----------------------------'''
##############################################################################################
debug('Stage: Processing restriction sites.')

## Process restriction site(s) as list
if args.restrictionEnzyme:
##    args.restrictionSite = rseqify_list([enzymes[e] for e in args.restrictionEnzyme.split(',')]) 
##    args.restrictionEnzyme = args.restrictionEnzyme.split(',')
    args.restrictionEnzyme, args.restrictionSite = rseqify_list2([e for e in args.restrictionEnzyme.split(',')])
    #enzymes = {args.restrictionEnzyme[i]:args.restrictionSite[i] for i in range(len(args.restrictionSite))}
    for i in range(len(args.restrictionSite)):
        enzymes[args.restrictionEnzyme[i]] = args.restrictionSite[i]
        enzymelookup[args.restrictionSite[i]].append(args.restrictionEnzyme[i])
    if len(args.restrictionEnzyme) <= 10:
        debug(str(zip(args.restrictionEnzyme, args.restrictionSite)))
else:
    args.restrictionSite = args.restrictionSite.split(',') 
    args.restrictionEnzyme = [e for e in args.restrictionSite]
    for e in args.restrictionEnzyme:
        enzymes[e]=e
        enzymelookup[e].append(e)





##############################################################################################
'''---------------------------[ PIPELINE FXN. ]----------------------------'''
##############################################################################################
def check():
    nseq=len(seqs)
    nhead = len(headers)
    debug('Loaded ' + str(nhead) + ' headers and ' + str(nseq) + ' associated sequences.')
    assert nseq == nhead

    


def pipeline(header, seq):
    ##############################################################################################
    '''---------------------------[ FIND RESTRICTION SITES ]----------------------------'''
    ##############################################################################################
    debug('Stage: Finding restriction sites.')

    debug('Identifying restriction site starts in sequence: ' + str(header) )
    if len(args.restrictionSite) > 1:
        debug('More than one enzyme detected - i.e. user provided more than one enzyme or recognition sequence.')
        
        # startsList is a list of lists
        startsList = [substringPositions1(seq, args.restrictionSite[i]) for i in range(len(args.restrictionSite))]

        n_startsList = len(startsList)
        n_ns_total = 0
        for i in range(n_startsList):
            n_sl_i = len(startsList[i])
            n_ns_total += n_sl_i
            debug('Detected ' + str(n_sl_i) + ' starts for ....' + args.restrictionSite[i])
        debug('Detected ' + str(n_ns_total) + ' sites overall.')
        debug('Non-empty starts lists came from:')
        if args.debug:
            for i in range(len(startsList)):
                if len(startsList[i])>0:
                    debug('\t'+args.restrictionSite[i])
                
        #bedlist = list of lists of tuples
        if args.bed:
            pass
            #bedlist = [(header, startsList[i][j], startsList[i][j]+len(args.restrictionSite[i]), args.restrictionSite[i]) for j in range(len(startsList)) for i in range(len(args.restrictionSite))]

        if args.digestTogether:
            if n_ns_total > 0: #Then combine list of lists into single list of starts and sort (as if only a single restrictionSite used)
                debug('DNA will be digested by more than one enzyme together in the same tube.')
                #starts = sorted( [x for e in startsList for x in e] )
                #starts = sorted( [x for i in range(len(startsList)) for x in startsList[i] ] )
                starts, starts_labels, sitelengths = zip(*sorted( [(x,args.restrictionEnzyme[i], len(args.restrictionSite[i])) for i in range(n_startsList) for x in startsList[i] ], key = lambda x: x[0] )) ## * means to unzip
                debug('Repeat: Detected ' + str(len(starts)) + ' sites overall.')
            else:
                starts, starts_labels, sitelengths = [], [], []


    elif len(args.restrictionSite) == 1:
        debug('Only one enzyme detected - i.e. user provided a single enzyme or recognition sequence.')

        # if only one restrictionSite, then just use a list of Starts
        starts = substringPositions1(seq, args.restrictionSite[0])
        starts_labels = [args.restrictionEnzyme[0]]*len(starts)
        sitelengths = [len(args.restrictionSite[0])]*len(starts)
        debug('Detected ' + str(len(starts)) + ' sites overall.')




    ##############################################################################################
    '''---------------------------[ CALCULATE LENGTHS ]----------------------------'''
    ##############################################################################################
    debug('Stage: Calculate lengths.')


    if args.digestTogether or len(args.restrictionSite) == 1: ## starts is a list
        debug('Single tube reaction.')
        if linear:
            debug('Calculating lengths given a linear molecule.')
            lengths, intervals = calculateLengths(starts, seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)
        elif args.circular:
            debug('Calculating lengths given a circular molecule.')
            lengths, intervals = calculateCircularLengths(starts, seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)

         
    else: ## starts is a list of lists
        ## Get lengths list for each list in list of lists
        ## Store as list of lists
        debug('More than one reaction.')
        lengthsList = []
        intervals = []
        for i in range(len(startsList)):
            starts = startsList[i]
            starts_labels = [args.restrictionEnzyme[i]]*len(starts)
            sitelengths = [len(args.restrictionSite[i])]*len(starts)
            if linear:
                debug('Calculating lengths given a linear molecule for reaction ' + str(i) + '.')
                lengths, beds = calculateLengths(startsList[i], seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)
                lengthsList += [lengths]
                intervals += [beds]
            elif args.circular:
                debug('Calculating lengths given a circular molecule for reaction ' + str(i) + '.')
                lengths, beds = calculateCircularLengths(startsList[i], seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)
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
        if len(args.restrictionSite) == 1 or args.digestTogether:
            for i in range(len(intervals)):
                interval = intervals[i]
                label = interval[2]
                labelseqs = ':'.join([enzymes[e] for e in label.split(':')])
                start_to_start_length = lengths[i]
                interval_length = interval[1]-interval[0]
                bed = [header, interval[0], interval[1], label+'_'+labelseqs, interval_length, '+', start_to_start_length]
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
            
        if args.rmap:
            print '\t'.join([str(e) for e in [header, len(seq), len(L)] + L])
        else:
            print ">"+header
            for i in range(len(L)):
                if args.proportion:
                    print args.mass*L[i]/float(len(seq))
                else:
                    print L[i]
    debug('Done with: '+header)
    debug('')



##############################################################################################
'''---------------------------[ Loading sequence(s) an running pipeline. ]----------------------------'''
##############################################################################################

## Get Sequence
debug('Stage: Loading sequence(s).')
if args.fasta:
    ## Traditional pipeline defined in restrictionDigest2.py
    debug('Processing sequence info from Fasta.')
    #headers, seqs = map(list, zip(*fastaToPython(args.fasta).iteritems()))
    for fa in SeqIO.parse(args.fasta, 'fasta'):
        pipeline(header=fa.description, seq=str(fa.seq))
    
elif args.commandline:
    debug('Processing sequence info from command-line.')
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
    pipeline(header, seq)






debug('All done.')










##############################################################################################
'''---------------------------[  ]----------------------------'''
##############################################################################################
##############################################################################################
'''---------------------------[  ]----------------------------'''
##############################################################################################
##############################################################################################
'''---------------------------[  ]----------------------------'''
##############################################################################################
##############################################################################################
'''---------------------------[  ]----------------------------'''
##############################################################################################




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


