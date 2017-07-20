#!/usr/bin/env python
###import xml.etree.cElementTree as ET
import lxml.etree
import lxml.builder    
import argparse

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Takes in options and writes the corresponding XML Protocol file for PBJelly.
    
EXAMPLE OUTPUT:
<jellyProtocol>
    <reference>/gpfs/scratch/jurban/pb/pbj/withallpbreads/pbout/hgap3.fasta</reference>  
    <outputDir>/users/jurban/scratch/pb/pbj/withpbreads/pbout</outputDir>
    <blasr>-minMatch 8 -sdpTupleSize 8 -minPctIdentity 65 -bestn 1 -nCandidates 10 -maxScore -500 -nproc 48 -noSplitSubreads</blasr>
    <input baseDir="/gpfs/scratch/jurban/pac_bio_data/filt/">
        <job>all_subreads.fastq</job>
    </input>
</jellyProtocol>

    
    """, formatter_class= argparse.RawTextHelpFormatter)

# required args
parser.add_argument('-r', '--reference', type=str, default=False, required=True, help='''Absolute path to reference fasta.''')
parser.add_argument('-o', '--outdir', type=str, default=False, required=True, help='''Absolute path to output dir.''')
parser.add_argument('-b', '--baseDir', type=str, default=False, required=True, help='''Absolute path to base directory where reads are stored.''')
parser.add_argument('-f', '--fastq', type=str, default=False, required=True, help='''Filename of fastq file to be used (expected to be in baseDir).''')
## rest have defaults
parser.add_argument('-m', '--minMatch', type=str, default='8', help='''blasr minMatch. Default: 8.''')
parser.add_argument('-s', '--sdpTupleSize', type=str, default='8', help='''blasr sdpTupleSize. Default: 8.''')
parser.add_argument('-p', '--minPctIdentity', type=str, default='75', help='''blasr minPctIdentity. Default: 75.''')
parser.add_argument('-B', '--bestn', type=str, default='1', help='''blasr bestn. Default: 1.''')
parser.add_argument('-N', '--nCandidates', type=str, default='10', help='''blasr nCandidates. Default: 10.''')
parser.add_argument('-M', '--maxScore', type=str, default='-500', help='''blasr maxScore. Default: -500.''')
parser.add_argument('-P', '--nproc', type=str, default='48', help='''blasr nproc. Default: 48.''')

args = parser.parse_args()



E = lxml.builder.ElementMaker()
JELLY = E.jellyProtocol
REF = E.reference
OUT = E.outputDir
BLASR = E.blasr
IN = E.input
IN_JOB = E.job

blasrline = "-minMatch " + args.minMatch + " -sdpTupleSize " + args.sdpTupleSize + " -minPctIdentity " + args.minPctIdentity + " -bestn " + args.bestn + " -nCandidates " + args.nCandidates + " -maxScore " + args.maxScore + " -nproc " + args.nproc + " -noSplitSubreads"

the_doc = JELLY(
        REF(args.reference),
        OUT(args.outdir),
        BLASR(blasrline),
        IN( IN_JOB(args.fastq), 
        baseDir=args.baseDir)
        )   

print lxml.etree.tostring(the_doc, pretty_print=True)
