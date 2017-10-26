#!/usr/bin/env python2.7
import sys, argparse, re
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Take in FASTA file.

Use regular expression to find sub-sequences.

Return FASTA file with each sub-sequence masked with selected character (default: N).

Use case: masking out the BNG restriction sites inside gaps with Ns to potentially help PBJelly work.
--> As it is now, PBJelly treats each sub-section of a gap separated by the BNG recognition sequence,
     as its own gap... Is this problematic?

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--fasta', '-f',
                   type=str, required=True,
                   help='''Path to input fasta file''')

parser.add_argument('--pattern', '-p',
                   type=str, default='NNCACGAGNN',
                   help='''After breaking scaffolds on Ns, minimum contig size to report. Default=NNCACGAGNN''')

parser.add_argument('--mask', '-m',
                   type=str, default='N',
                   help='''Character used to mask sub-sequences. Default = N.''')

parser.add_argument('--coords', '-c',
                   action='store_true', default=False,
                   help='''At end, write BED file of coordinates of masked sequences.''')

args = parser.parse_args()



pattern = re.compile(args.pattern)


coords = {}
for fa in SeqIO.parse(args.fasta, 'fasta'):
    seq = str(fa.seq)
    newseq = ''
    idx = 0
    cnt = 0
    coords[fa.name] = []
    for e in re.finditer(pattern, seq):
        cnt += 1
	newseq += seq[idx:e.start()] + args.mask * (e.end()-e.start())
	idx = e.end()
	coords[fa.name] += [(e.start(), e.end())]
    newseq+=seq[idx:]
    print ">"+fa.name + "_" + str(cnt) + "masks"
    print newseq

if args.coords:
    outpre = (".").join( args.fasta.split("/")[-1].split(".")[:-1] )
    out = open(outpre+'.bed', 'w')
    for seq in sorted(coords.keys()):
        for pos in coords[seq]:
            out.write( ("\t").join( [str(e) for e in (seq, pos[0], pos[1])] ) + "\n" )

    
