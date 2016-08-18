#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

## DEFINE PARSER ARGUMENTS
parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in canu fasta(s) and outputs only those satisfying some filter criteria with info in canu header..
    Example Canu header:
    tig00000000 len=25777 reads=60 covStat=36.14 gappedBases=no class=contig suggestRepeat=no suggestCircular=no
    
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('canu', metavar='canu', nargs='+',
                   type= str, 
                   help='''Path(s) to canu fasta file(s). Can handle more than one, though it might not be recommended to use more than one.''')
parser.add_argument('-r', '--minreads', type=int, default=0,
                    help='''Only return fasta entries with this many reads or more.''')
parser.add_argument('-R', '--maxreads', type=int, default=1000000000,
                    help='''Only return fasta entries with this many reads or less.''')

args = parser.parse_args()



i=0
j=0
for canu in args.canu:
    i += 1
    for fa in SeqIO.parse(canu, "fasta"):
        j += 1
 
        desc = fa.description.split()
        length = float(desc[1].split("=")[1])
        nreads = float(desc[2].split("=")[1])
        covstat = float(desc[3].split("=")[1])
        gapped = desc[4].split("=")[1]
        seqclass = desc[5].split("=")[1]
        repeat = desc[6].split("=")[1]
        circular = desc[7].split("=")[1]
        if nreads >= args.minreads and nreads <= args.maxreads:
            print ">"+fa.description
            print str(fa.seq)
        
        
