#!/usr/bin/env python
import argparse, os
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
parser.add_argument('--start', '-s',
                   type=int, default=1, required=False,
                   help='''When files are split, they are given names modified by adding subsequent integers starting from 1 by default.
Use this flag to provide an alternative integer to start with -- e.g. 0 if you need/want it to be 0-based.''')

parser.add_argument('--outdir', '-o', type=str, default=False,
                   help='''Provide path to place split fasta files. Default: same directory original file is in.''')

parser.add_argument('--wd', '-w', action='store_true', default=False,
                   help='''This is a shortcut for -o to place split fasta files in working directory (dir script is used in). Same as "-o ./"''')


args = parser.parse_args()

if args.wd:
    args.outdir = "./"

## functions
def get_base_name(fasta_handle, outdir=False):
    fa_path = os.path.dirname(fasta_handle)
    fa_base = (".")join( os.path.basename(fasta_handle).split(".")[:-1] )
    if outdir:
        if outdir.endswith("/"):
            base_name = outdir + fa_base
        else:
            base_name = outdir + "/" + fa_base
    else:
        base_name = fa_path + "/" + fa_base
    return base_name


def checkdir(outdir):
    if not os.path.exists(outdir):
        os.system( "mkdir " + outdir )

def split(fasta_handle, numentries, base_name, file_num = 1):
    count = 0
    filename = base_name + "." + str(file_num) + ".fa"
    fh = open(filename, 'w')
    for fa in SeqIO.parse(fasta_handle, "fasta"):
        if count == numentries:
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
base_name = get_base_name(args.fasta, args.outdir)
checkdir(args.outdir)
split(args.fasta, args.numentries, base_name, args.start)
        
