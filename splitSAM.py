#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
import pysam

parser = argparse.ArgumentParser(description="""

DESCRIPTION - version 21Nov2016
    Takes in SAM or BAM.
    For now only outputs BAM.
    Default split is into 8 files.
    For now it is required that you provide number of reads/alignments in the input SAM/BAM.
        You can do this with 'samtools view -c in.bam'.
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--sam', '-s',
                   type= str, default=False,
                   help='''Input file in SAM format.''')
parser_input.add_argument('--bam', '-b',
                   type= str, default=False,
                   help='''Input file in BAM format.''')
parser.add_argument('--nfiles',type=int, default=10,
                    help='''Split into this many files. Provide Int. Default: 8.''')
parser.add_argument('--nreads',type=int, required=True,
                    help='''For now, provide number of reads/alignments in input SAM/BAM file.
                            Do this by 'samtools view -c file.bam'.
                            This is a required argument and there is no default.''')
args = parser.parse_args()



if args.sam:
    samfile = pysam.Samfile(args.sam, "r")
    pre = args.sam.strip().split('.')[:-1]
elif args.bam:
    samfile = pysam.Samfile(args.bam, "rb")
    pre = args.bam.strip().split('.')[:-1]



nreadsperfile = (args.nreads//args.nfiles)+1
nreadsinfile = 0

#initialize
filenum=1
outsam = pysam.AlignmentFile(pre+ str(filenum)+".bam", "wb", template=samfile)
# begin iter loop
for read in samfile.fetch():
    outsam.write(read)
    nreadsinfile += 1
    if nreadsinfile == nreadsperfile:
        outsam.close()
        filenum += 1
        nreadsinfile = 0
        outsam = pysam.AlignmentFile(pre+ str(filenum)+".bam", "wb", template=samfile)
#close last open file
outsam.close()



