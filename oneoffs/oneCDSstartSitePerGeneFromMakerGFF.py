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



Goal:
Get single CDS start site entry for each gene.
For genes on + strand, this will be lowest integer in CDS coords.
For genes on - strand, this will be the highest integer in CDS coords.
It will be reported as BED.



Example + strand:
Input:
contig_1        maker   CDS     25763   25788   .       +       0       ID=Bcop_v1_g022643-RA:cds;Parent=Bcop_v1_g022643-RA
contig_1        maker   CDS     25832   26014   .       +       1       ID=Bcop_v1_g022643-RA:cds;Parent=Bcop_v1_g022643-RA
contig_1        maker   CDS     26140   26272   .       +       1       ID=Bcop_v1_g022643-RA:cds;Parent=Bcop_v1_g022643-RA
contig_1        maker   CDS     28452   28526   .       +       0       ID=Bcop_v1_g022643-RA:cds;Parent=Bcop_v1_g022643-RA

Output:
contig_1        maker   cdsStart     25763   25763   .       +       0       ID=Bcop_v1_g022643.....etc
OR
contig_1    25762    25763    Bcop_v1_g022643    .    +

Example from - strand:
Input:
contig_4        maker   CDS     2701    2853    .       -       0       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     3471    3613    .       -       2       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     10455   10596   .       -       0       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     10692   10892   .       -       0       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     10954   11036   .       -       2       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     13261   13330   .       -       0       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     13732   13902   .       -       0       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     13968   14109   .       -       1       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA
contig_4        maker   CDS     20942   21219   .       -       0       ID=Bcop_v1_g007771-RA:cds;Parent=Bcop_v1_g007771-RA

Output:
contig_4        maker   cdsStart     21219   21219   .       -       0       ID=Bcop_v1_g007771...etc
OR
contig_4    21218    21219    Bcop_v1_g007771    .    -

    """, formatter_class= argparse.RawTextHelpFormatter)




parser.add_argument('gff', metavar='gff', nargs=1,
                   type= str, 
                   help='''Path to Maker2 GFF.''')


parser.add_argument('-B', '--bedout',
                   action='store_true', default=False,
                   help='''
                        Return BED entries rather than GFF entries.
                        
                        ''')
args = parser.parse_args()





def read_gff(fh):
    gff = GFF()
    with open(fh) as gff_fh:
        for line in gff_fh:
            gff.add(line)
    return gff



def run(args):
    gff = read_gff(args.gff[0])
    for gff_entry in gff.get_single_CDS_start_site_per_gene():
        if args.bedout:
            gff_entry.print_as_bed()
        print gff_entry






run(args)
        

