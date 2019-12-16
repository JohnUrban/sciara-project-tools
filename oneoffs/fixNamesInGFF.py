#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - fix values in GFF description given a map/dictionary file.

Input:
    - GFF (e.g. from Maker2)
    - map (e.g. from maker_map_ids)
        - two col file: old_name new_name

Output:
    - Updated GFF where old names are converted to new names

Use case:
    - I found that when trying to check my GFF with table2asn_GFF from NCBI
        - it gave errors like:
            - Bad data line: Record references non-existant Parent=maker-contig_38-snap-gene-0.419-mRNA-2
        - upon closer inspection these were coming up in the Parent= part of the description.
            - contig_39	maker	exon	1896096	1896162	.	-	.	ID=Bcop_v1_g017328-RC:exon:1305;Parent=Bcop_v1_g017328-RC,maker-contig_39-snap-gene-0.919-mRNA-3
        - So I guess the maker_map_ids script did not work for this...



STOPPED DEVELOPMENT.
- I found another source to the issue


    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--gff', '-g', '-i', '-f',
                   type= str, default=False, required=True,
                   help='''Path to GFF.''')

parser.add_argument('--map', '-m',
                   type= str, default=False, required=True,
                   help='''Path to map file.''')

args = parser.parse_args()




# Build a class for the dictionary that returns values not already present
class MapDict(object):
    def __init__(self, d):
        self.d = d
    def get_value(self, key):
        try:
            self.d[key]
        except KeyError:
            return key
        
# Build dictionary from map file
d = {}
with open(args.map) as f:
    # Read in all lines
    lines = f.readlines()

    # Go through lines.
    for line in lines:
        # Process line
        line = line.strip().split('\t')

        # Add oldnames (col1) as key and newnames (col2) as value.
        d[line[0]] = line[1]

# Then MapDict object from d
mapd = MapDict(d)





# Process GFF
with open(args.gff) as f:
    for line in f:
       pass

