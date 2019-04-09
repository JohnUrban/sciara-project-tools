#!/usr/bin/env python2.7
import sys
import argparse
import random
import os
#sys.path.append("/Users/johnurban/searchPaths/pyScripts/")
from pfmTools import *

VersionInfo="""Version 1.0
This version was created on 10/23/2013 by John Urban."""

parser = argparse.ArgumentParser(description=VersionInfo+"""

    DESCRIPTION

        Takes in fasta file output from BEDtools fastaFromBed.
        All records need to be the same length.
        Returns a file with a position frequency matrix.

        The output file has a single header line where each column is labeled the
        nucleotide symbol counted in that column. Each line/row in the file is a
        position along the sequences from 1:N.


    EXAMPLE OUTPUT:
    
        A\tC\tG\tT
        100\t20\t10\t100
        .\t.\t.\t.
        .\t.\t.\t.
        .\t.\t.\t.
        10\t100\t100\t20

    NOTES:
        When piping into this, it first writes out the stdin as temporary file (that is removes at the end).
        What is piped in should be a fasta file -- most likely the output of fastaFromBed.

        Creating fasta:
        The following starts with a 1bp summits BED file, extends it, keeps only those that did not fall of chr and do not overlap gaps, and makes fasta.
        slopBed -i summits.bed -g genome -b size | awk '{if ($3-$2 == X) print $0}' | intersectBed -v -a - -b gaps | fastaFromBed -fi fasta -bed - -fo out.fasta
 
        Fasta file format is completely flexible.
            i.e. it does not have to be that a sequence is all on 1 line
            i.e. it does not have to be the case that all the formatting (line widths for each fasta sequence entry) are uniform
        Output from BEDtools fastaFromBed will be both on 1 line and uniform though.
        If you feel you must format your fasta file, use fasta_formatter from fastX toolkit.
            


    TODO:
        -- Make it possible to use -o
            -- Right now it only does std out
            -- For most things this should be fine. Just redirect into a file.
        -- I think this is actually slower than the matlab version.
            -- not a direct comparison though, b/c I ran matlab on Oscar and this on my MBP
            -- A major difference is that the matlab script loads it all into memory which is possible on Oscar
            -- The python script goes line by line assembling and analyzing 1 sequence at a time
            -- Thus its memory footprint is smaller, but is a little slower
                -- I dont think I was able to ever actually use the matlab version on my personal computer
            -- Benchmarks:
                196705 sequences, each 4001 bp long, ~787 Mbp total, ~6 m 30 s
                161933 sequences, each 4001 bp long, ~647.9 Mb total, 5 m 20 s
                66750 sequences, each 4001 bp long, ~267.1 Mb total, ~2 m 5 s
        -- Have it also write out the number of sequences analyzed
            -- This will be used in R
            -- Alternatively, have the option to do Pfm, Pwm, or both
                    -- where pwm = pfm/numSeq
    

    """,
    formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-fi', '--fasta',
                    type=str,
                    help='''The fasta file to analyze. 
This is often the output of BEDtools fastaFromBed.
All fasta entries should be the same length.
Can also pipe in the output -- use "-" or "stdin" instead.
                    ''',
                    required=True)

parser.add_argument('-o', '--outName',
                    type=str,
                    help='''The name of the output file (optional).
Default is stdout.
                    ''',
                    default="stdout")
parser.add_argument('-c', '--consensus',
                    help="Return consensus instead of PFM. Consensus just finds the symbol with max counts for each position. ",
                    action="store_true")
parser.add_argument('-v', '--verbose',
                    help="Use to increase verbosity while program is running",
                    action="store_true")
parser.add_argument('--version', action='version', version=VersionInfo)

args = parser.parse_args()


#"---------------------------[ Check and pare arguments ]----------------------------"
stdin=False
if args.fasta == '-' or args.fasta == 'stdin':
    stdin=True
    outSeed = random.randint(1000000,9999999)
    args.fasta = "tempFasta." + str(outSeed) + ".txt"
    fastaOut = open(args.fasta,'w')
    line = sys.stdin.readline()
    while line != '':
        fastaOut.writelines(line)
        line = sys.stdin.readline()
    fastaOut.close()

#"---------------------------[ Functions ]----------------------------"



#'-----------------[ execution ]----------------------------------'
pfmDict = pfmDictFromFasta(args.fasta)
if not args.consensus:
    pfmFromPfmDict(pfmDict)
else:
    pfmConsensus(pfmDict, args.verbose)

            
#'-----------------[ Clean up ]----------------------------------'
if stdin:
    os.remove(args.fasta)
