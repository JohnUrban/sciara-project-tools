#!/usr/bin/env python2.7
import sys, re, string
import argparse
import random
import os
from itertools import product
from Bio import SeqIO
from collections import defaultdict
from restrictionDigestTools import *
from restrictionEnzymes import *

VersionInfo="""Version beta 3.1

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

parser.add_argument('-re', '--regexpipeline',
                    action='store_true',
                    default=False,
                    help='''Use version 3's regex pipeline. Might be faster. Definitely many times faster when using enzymes with N/R/Y in sequence.
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



verbose=(args.debug or args.verbose)

##########################################################################
'''---------------------------[ Messages that no longer precede their opertions ]----------------------------'''
##########################################################################

debug('Stage: Defining functions.',args)
debug('Defining enzymes and restriction sites.',args)

        


if args.listEnzymes:
    print
    for key in sorted(enzymes.keys()):
        regex = ''
        if enzymes[key] != '*':
            regex = regexify(enzymes[key])
        print key + "\t" + enzymes[key] + "\t" + regex
    print
    quit()





##############################################################################################
'''---------------------------[ Loading sequence(s) an running. ]----------------------------'''
##############################################################################################


## Get Sequence
debug('Stage: Loading sequence(s).',args)
args.sites_processed = False ## This needs to be set as False here to trigger restriction site processing, which then sets it to false to avoid re-doing for all subsequent fasta entries

if args.fasta:
    debug('Processing sequence info from Fasta.',args)
    if args.regexpipeline:
        PPLN=regex_ver2pipeline
    else:
        ## Traditional pipeline defined in restrictionDigest2.py
        #headers, seqs = map(list, zip(*fastaToPython(args.fasta).iteritems()))
        PPLN=ver2pipeline
    for fa in SeqIO.parse(args.fasta, 'fasta'):
        PPLN(header=fa.description, seq=str(fa.seq), args=args)
    
elif args.commandline:
    debug('Processing sequence info from command-line.',args)
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
    debug('Loaded ' + str(header) + ': ' + str(seqlen) + ' bp.',args)
    pipeline(header, seq,args)






debug('All done.',args)

