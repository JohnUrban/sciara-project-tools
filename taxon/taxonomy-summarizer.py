#!/usr/bin/env python

##import sys
import argparse
import gzip
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

DESCRIPTION
  Takes in tab-separated value (tsv) file output of cnv-taxonomyFromTaxID.py:
  (1) name of sequence being characterized
  (2) Closest level to Insecta
  (3) species info

  Outputs 4 files.
  (1) Insecta_hits
      Contains seq names that had at least Xi Insecta characterizations (default: Xi = 1)
  (2) Metazoa_hits
      Contains seq names that had xi < Xi, but at least Xm Metazoa characterizations (default: Xm = 1)
  (3) Eukaryota_hits
      Contains seq names that had xi < Xi and xm < Xm, but at least Xe Eukaryota characterizations (default: Xe = 1)
  (4) Other_hits
      Contains seq names that had xi < Xi, xm < Xm, and xe < Xe.

  All output files contain seqname in column1 followed by:
  column 2 = Insecta:xi
  column 3 = Metazoa:xm
  column 4 = Eukaryota:xe
  column 5 = Archaea:xa
  Column 6 = Bacteria:xb
  Column 7 = Other:xo
  Columns 8 to N = all colon-separated species:count found, "-" if None.
    """, formatter_class= argparse.RawTextHelpFormatter)
##  Column 7 = comma-separated list of all colon-separated insecta_species:count found, "-" if None.
##  Column 8 = comma-separated list of all colon-separated Metazoa_species:count found, "-" if None.
##  Column 9 = comma-separated list of all colon-separated Eukaryota_species:count found, "-" if None.
##  Column 10 = comma-separated list of all colon-separated Archaea_species:count found, "-" if None.
##  Column 11 = comma-separated list of all colon-separated Bacteria_species:count found, "-" if None.
##  Column 12 = comma-separated list of all colon-separated other_things_found:count, "-" if None.
parser.add_argument('--inputfile', '-i',
                   type=str, default=False, required=True,
                   help='''Provide tab-separated file output of cnv-taxonomyFromTaxID.py''')
parser.add_argument('--outprefix', '-o',
                   type=str, default=False, required=True,
                   help='''Provide output prefix.''')
parser.add_argument('--xi', '-xi',
                   type=int, default=1,
                   help='''Provide minimum number of Insecta hits to be included in Insecta file. Default = 1.''')
parser.add_argument('--xm', '-xm',
                   type=int, default=1,
                   help='''Provide minimum number of Metazoa hits to be included in Metazoa file (if does not meet reqs for Insecta). Default = 1.''')
parser.add_argument('--xe', '-xe',
                   type=int, default=1,
                   help='''Provide minimum number of Eukaryota hits to be included in Eukaryota file (if does not meet Insecta nor Metazoa reqs). Default = 1.''')
args = parser.parse_args()


########################################################################################################
########################################################################################################
########################################################################################################
## FUNCTIONS
########################################################################################################
########################################################################################################
########################################################################################################


def open_file_connection(fh):
    ## Open Connection to gz or non-gz file
    if fh == "-" or fh == "stdin":
        connection = sys.stdin
    else:
        connection = gzip.open(fh)
        try:
            connection.next()
            connection.seek(0)
        except IOError:
            connection = open(fh)
    return connection


def get_line(fh):
    try: 
        line = fh.next()
        line = line.strip().split("\t")
        return line
    except StopIteration:
        return False

    
########################################################################################################
########################################################################################################
########################################################################################################
## EXECUTE
########################################################################################################
########################################################################################################
########################################################################################################


fh = open_file_connection(args.inputfile)

proximity2insecta = defaultdict(lambda: defaultdict(int))
species_counts = defaultdict(lambda: defaultdict(int))

for line in fh:
    line = line.strip().split()
    proximity2insecta[line[0]]
    proximity2insecta[line[0]][line[1]] += 1
    species_counts[line[0]]
    species_counts[line[0]][line[2]] += 1
fh.close()

fi = open(args.outprefix+".insecta.out", 'w')
fm = open(args.outprefix+".metazoa.out", 'w')
fe = open(args.outprefix+".eukaryota.out", 'w')
fo = open(args.outprefix+".other.out", 'w')

for seqname in proximity2insecta.keys():
    ## keys will be seq names
    total = sum(proximity2insecta[seqname].values())
    xi = proximity2insecta[seqname]["Insecta"]
    xm = proximity2insecta[seqname]["Metazoa"]
    xe = proximity2insecta[seqname]["Eukaryota"]
    xa = proximity2insecta[seqname]["Archaea"] #will enter this name and 0 into dict if not present
    xb = proximity2insecta[seqname]["Bacteria"] #will enter this name and 0 into dict if not present
    proximity2insecta[seqname]["Other"] = total - xi-xm-xe-xa-xb ## Other is a redundant count -- the other thing(s) will be counted by w/e is in that field in above loop
    outline = ("\t").join([seqname] + [str(k)+":"+str(v) for k,v in proximity2insecta[seqname].iteritems()] + [str(k)+":"+str(v) for k,v in species_counts[seqname].iteritems()])
    if xi >= args.xi:
        fi.write(outline+"\n")
    elif xm >= args.xm:
        fm.write(outline+"\n")
    elif xe >= args.xe:
        fe.write(outline+"\n")
    else:
        fo.write(outline+"\n")

            
fi.close()
fm.close()
fe.close()
fo.close()

