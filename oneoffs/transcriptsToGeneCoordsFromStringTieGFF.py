#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict
import numpy as np
from gfftools import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION -
Version 0.1-20190529

GFF3: https://useast.ensembl.org/info/website/upload/gff3.html
seqid - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
source - name of the program that generated this feature, or the data source (database or project name)
type - type of feature. Must be a term or accession from the SOFA sequence ontology
start - Start position of the feature, with sequence numbering starting at 1.
end - End position of the feature, with sequence numbering starting at 1.
score - A floating point value.
strand - defined as + (forward) or - (reverse).
phase - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
attributes - A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation for more details.


StringTie GFF3 Example for 3 transcripts from 1 gene:
contig_103      StringTie       transcript      455301  456207  1000.00 -       .       ID=MSTRG.41.1;geneID=MSTRG.41
contig_103      StringTie       exon    455301  455938  1000.00 -       .       Parent=MSTRG.41.1
contig_103      StringTie       exon    455994  456207  1000.00 -       .       Parent=MSTRG.41.1
contig_103      StringTie       transcript      455301  456207  1000.00 -       .       ID=MSTRG.41.2;geneID=MSTRG.41
contig_103      StringTie       exon    455301  455930  1000.00 -       .       Parent=MSTRG.41.2
contig_103      StringTie       exon    455994  456207  1000.00 -       .       Parent=MSTRG.41.2
contig_103      StringTie       transcript      455316  456191  1000.00 -       .       ID=MSTRG.41.3;geneID=MSTRG.41
contig_103      StringTie       exon    455316  455938  1000.00 -       .       Parent=MSTRG.41.3
contig_103      StringTie       exon    456002  456191  1000.00 -       .       Parent=MSTRG.41.3


Goal:
Get single entry for gene that marks the longest span.

Here:
contig_103      StringTie       gene      455301  456207  . -       .       ID=MSTRG.41;geneID=MSTRG.41



    """, formatter_class= argparse.RawTextHelpFormatter)




parser.add_argument('gff', metavar='gff', nargs=1,
                   type= str, 
                   help='''Path to StringTie GFF.''')

args = parser.parse_args()





def read_gff(fh):
    gff = StringTieGFF()
    with open(fh) as gff_fh:
        for line in gff_fh:
            gff.add(line)
    return gff



def run(fh):
    gff = read_gff(fh)
    for line in gff.get_merged_gene_entries():
        print line






run(args.gff[0])
        
