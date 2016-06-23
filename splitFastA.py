#!/usr/bin/env python
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Take a multi-fasta file (.fa file with >1 entry) and split it into multiple subfiles
    by specifying how many entries per file. The last file to be made will often have less than this number.
    
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--fasta', '-f', required=True,
                   type=str, default=False,
                   help='''Provide path to fasta file''')
parser.add_argument('--numentries', '-n',
                   type=int, default=False, required=True,
                   help='''Provide number of entries to put in each split file.''')

args = parser.parse_args()


## functions
def split(fasta_handle):
    base_name = fasta_handle.split(".")[0]
    file_num = 1
    count = 0
    filename = base_name + "." + str(file_num) + ".fa"
    fh = open(filename, 'w')
    for fa in SeqIO.parse(fasta_handle, "fasta"):
        if count == args.numentries:
            fh.close()
            count = 0
            file_num += 1
            filename = base_name + "." + str(file_num) + ".fa"
            fh = open(filename, 'w')
        count += 1
        fh.write(">"+str(fa.id)+"\n")
        fh.write(str(fa.seq)+"\n")
    fh.close()



## Execute
split(args.fasta)
        
