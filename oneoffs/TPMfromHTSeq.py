#!/usr/bin/env python2.7
import sys, os, argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Given a transcript length file and any number of transcript count files, convert counts to TPM.

1. Transcript length file is tab-sep name and length.
2. Transcript count file is tab-sep name and count.

Count can be from HTSeq-count.
With dU strand-specific method where R1-R2 give the RF pattern I've done this command.
hisat2 -p ${P} -x ${HIDX} -1 $R1 -2 $R2 --rna-strandness RF | samtools sort > ${PRE}.bam
htseq-count -f bam -r pos -s reverse -t exon -i ID ${BAM} Bradysia_coprophila.Bcop_v1.0_gene_set.gff 1> htseq-count.${PRE}.txt 2>htseq-count.${PRE}.err 

That gives Exon level counts.
I am using a Maker2 GFF3.
Thus Exon level names look like:
Bcop_v1_g000001-RA:exon:3042

Transcript level looks like:
Bcop_v1_g000001-RA

Gene level looks like:
Bcop_v1_g000001

To get transcript level counts from just exons, do:
awk '{gsub(/:/,"\t"); a[$1]+=$4}END{for (e in a) print e"\t"a[e]}' $f > ${PRE}.transcript-level.txt

To get gene level counts from just exons, do:
awk '{gsub(/:/,"\t"); sub(/-/,"\t"); a[$1]+=$5}END{for (e in a) print e"\t"a[e]}' $f > ${PRE}.gene-level.txt


To go from exons to transcript level TPMs, simply do the above to get transcript counts, then submit here to get transcript TPMs.
To go from exons to gene-level TPMs, get transcript-level TPMs then sum those for each gene with awk similar to above.
Alternatively, pick a single length (e.g. longest transcript) for each gene.


Because HTSeq has non-genes such as __not_aligned, __no_feature, etc -- you might want to remove those beforehand.

Try:
grep -v "^__" file | TPMfromHTSeq.py -l lengths.txt -
or
awk '$1 !~ /^__/' file | TPMfromHTSeq.py -l lengths.txt -

""", formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('counts', metavar='counts', nargs='*',
                   type= str, 
                   help='''Paths to as many tables as you want. They will be treated as one giant table.
                        This can also be left empty for standard in or specified as - or stdin.''')

parser.add_argument('-l', '--lengths', type=str, required=True,
                    help='''Path to file of transcript lengths.''')


args = parser.parse_args()




## Functions
def table_to_dict(fh):
    return dict((e[0], int(e[1])) for e in [line.strip().split() for line in fh.readlines()])

def read_table(f):
    if f in ('stdin', '-') or f.startswith('<('):
        fh = sys.stdin
        return table_to_dict(fh)
    else:
        with open(f) as fh:
            return table_to_dict(fh)


def TPM(counts, lengths, valtype=str):
    norm = dict((name, counts[name]/float(lengths[name])) for name in counts.keys())
    scale = sum(norm.values())
    return dict((name, valtype(1e6*norm[name]/scale)) for name in norm.keys())

def process_counts(countsfilepaths):
    # Read in lengths
    lengths = read_table(args.lengths)
    
    for f in countsfilepaths:
        # Read in count
        counts = read_table(f)

        # Create output
        print '\n'.join('\t'.join(e) for e in sorted(TPM(counts,lengths).iteritems()))


def main(args):
    # Loop over count files
    process_counts(args.counts)




## Execute
try:
    main(args)
except IOError:
    pass
                

