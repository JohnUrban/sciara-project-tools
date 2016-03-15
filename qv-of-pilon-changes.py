#!/usr/bin/env python
import gzip, sys, argparse
import numpy as np
from Bio import SeqIO
from assembly_class import *
from pilon_change_class import *
from helper_functions import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    This script assumes one had a genome assembly that was polished by quiver,
    and that quiver assembly was further polished by Pilon.

    It takes the pilon.changes file and the fastq file of the assembly output by quiver.

    For each change, it reports the median quality value of the bases the were changed.

    The output is the original file with additional columns.

    Note to self:
    All Pilon genome coordinates are 1-based. Python is 0-based.
    Do proper accounting.
    Same as converting BED to SAM.
    
    TODO:
    In future - perhaps also reports coverage of pacbio reads over the bases.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-p", "--pilonfile",
                   type= str, default=False, required=True,
                   help='''Input pilon.changes file.''')
parser.add_argument("-q", "--quiverfile",
                   type= str, default=False, required=True,
                   help='''Input quiver-polished-assembly.fastq file.''')
parser.add_argument("-v", "--verbose",
                   action = "store_true", default=False,
                   help='''Talk about what's going on.''')
args = parser.parse_args()






##execute
changes = open(args.pilonfile, "r")

if args.verbose:
    sys.stderr.write("Loading assembly....\n")
assembly = Assembly(args.quiverfile,fastx="fastq")
assembly.load_assembly()

if args.verbose:
    sys.stderr.write("Comparing changes to assembly...\n")

characterize_changes(changes, assembly)












### TESTS
def does_seq_from_assembly_match(changes):
    for line in changes:
        change = Change(line)
        seq = change.get_old_seq()
        contig, start, end = change.get_coord()
        assemb_seq = assembly.extract_info(contig,start,end)
        assemb_seq = assemb_seq.seq
        assemb_seq = str(assemb_seq)
        if seq != ".":
            if not assembly.is_correct_seq(contig, start, end, seq, makeupper=True):
                print seq
                print assemb_seq
                print

