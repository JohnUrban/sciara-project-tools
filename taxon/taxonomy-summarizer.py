#!/usr/bin/env python

import sys
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

  All output files contain seqname in column1 followed by taxon levels and numbers of hits....
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
parser.add_argument('--allseqnames', '-a', default=False,
                    help=''' Optional. Provide file with all seq names you are attempting to classify, 1 per line (e.g. all contig names in assembly).
                    In a blast hits file (and therefore taxonomy output file use as input here), not all sequences you blasted necessarily have hits.
                    Therefore, by default they will not be classified as anything and would not be part of the output here.
                    This option just allows you to report any sequence name not in given list that was not classified -- reported in 'other' file.
                    ''')
parser.add_argument('--outprefix', '-o',
                   type=str, default=False, required=True,
                   help='''Provide output prefix.''')
##parser.add_argument('--xi', '-xi',
##                   type=int, default=1,
##                   help='''Provide minimum number of Insecta hits to be included in Insecta file. Default = 1.''')
##parser.add_argument('--xm', '-xm',
##                   type=int, default=1,
##                   help='''Provide minimum number of Metazoa hits to be included in Metazoa file (if does not meet reqs for Insecta). Default = 1.''')
##parser.add_argument('--xe', '-xe',
##                   type=int, default=1,
##                   help='''Provide minimum number of Eukaryota hits to be included in Eukaryota file (if does not meet Insecta nor Metazoa reqs). Default = 1.''')

parser.add_argument('--superkingdom', type=str, default="Eukaryota",
                    help='''Provide super kingdom of your target species.
Default: Eukaryota.
Example of other: Bacteria, Archaea.''')

parser.add_argument('--kingdom', type=str, default="Metazoa",
                    help='''Provide kingdom of your target species. Default: Metazoa.''')


parser.add_argument('--phylum', type=str, default="Arthropoda",
                    help='''Provide phylum of your target species. Default: Eukaryota.''')

parser.add_argument('--Class', type=str, default="Insecta",
                    help='''Provide class of your target species. Default: Insecta.''')


parser.add_argument('--order', type=str, default="Diptera",
                    help='''Provide order of your target species. Default: Diptera.''')


parser.add_argument('--family', type=str, default="Sciaridae",
                    help='''Provide order of your target species. Default: Sciaridae.''')

parser.add_argument('--genus', type=str, default="Bradysia",
                    help='''Provide order of your target species. Default: Bradysia.''')

parser.add_argument('--species', type=str, default="Bradysia coprophila",
                    help='''Provide order of your target species. Default: Bradysia coprophila.''')

parser.add_argument('--levels', type=str, default='species,genus,family,order,class,phylum,kingdom,superkingdom,othersuperkingdoms',
                    help='''Use only the given levels for classification. Give comma-separated list. Default (all): species,genus,family,order,class,phylum,kingdom,superkingdom,othersuperkingdoms''')

parser.add_argument('--level_mins', type=str, default='1,1,1,1,1,1,1,1,1',
                    help='''Using only the given levels for classification, assign a sequence to a file if it has at least given number of hits.
Give comma-separated list - the same length as --levels. Default (all): 1,1,1,1,1,1,1,1,1''')

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
    line = line.strip().split("\t")
    proximity2insecta[line[0]]
    proximity2insecta[line[0]][line[1]] += 1
    species_counts[line[0]]
    species_counts[line[0]][line[2]] += 1
fh.close()

if args.allseqnames:
    with open(args.allseqnames, 'r') as fh:
        for line in fh:
            proximity2insecta[line.strip().split()[0]] ## allows files with more than 1 col, but names must be first

all_levels = 'species,genus,family,order,class,phylum,kingdom,superkingdom'.split(',')
all_level_args = [args.species, args.genus, args.family, args.order, args.Class, args.phylum, args.kingdom, args.superkingdom]
levels = []
level_vals = []
mins = [int(e) for e in args.level_mins.split(",")]
level_mins = []
for i in range(len(all_levels)):
    if all_levels[i] in args.levels:
        levels.append(all_levels[i])
        level_vals.append(all_level_args[i])
        level_mins.append(mins[i])
if "othersuperkingdoms" in args.levels:
    other_super_kingdoms = ["Eukaryota", "Bacteria", "Archaea"]
    other_super_kingdoms.remove(args.superkingdom)
    for i in range(len(other_super_kingdoms)):
        levels.append("othersuperkingdom"+str(i+1))
        level_vals.append(other_super_kingdoms[i])
        level_mins.append(mins[-1])

files = {}
fnames = {}
for i in range(len(levels)):
    fname = args.outprefix + "." + levels[i] + "." + ("_").join(level_vals[i].split()) + ".out"
    fnames[levels[i]] = fname
    files[levels[i]] = open(fname, 'w')
fname = args.outprefix + ".other.out"
other = open(fname, 'w')
fnames["Other"] = fname
##fname = args.outprefix + ".unclassified.out"
##other = open(fname, 'w')
##fnames["unclassified"] = fname


for seqname in proximity2insecta.keys():## keys are seq names
    ## get total number of hits
    total = sum(proximity2insecta[seqname].values())
    proximity2insecta[seqname]["Other"] = total

    ## go through levels to account for number of hits not accounted for in levels
    for i in range(len(levels)):
        proximity2insecta[seqname]["Other"] -= proximity2insecta[seqname][level_vals[i]]

    ## Construct outline
    outline = ("\t").join([seqname] + [str(k)+":"+str(v) for k,v in proximity2insecta[seqname].iteritems() if v > 0] + [str(k)+":"+str(v) for k,v in species_counts[seqname].iteritems() if v > 0]) + "\n"

    # go through levels, write info into first level file that meets rules for proximity to sp of interest (level mins)
    level_found=False
    for i in range(len(levels)):
        if proximity2insecta[seqname][level_vals[i]] >= level_mins[i]:
            level_found=True
            files[levels[i]].write(outline)
            break
        
    if not level_found:
        other.write(outline)
            
## close out files
import os
for i in range(len(levels)):
    files[levels[i]].close()
other.close()

# remove empty files
for fname in fnames.keys():
    if os.stat(fnames[fname]).st_size == 0:
        sys.stderr.write("Nothing written to " + fnames[fname] + "... removing...\n")
        os.remove(fnames[fname])
##
##fi = open(args.outprefix+".insecta.out", 'w')
##fm = open(args.outprefix+".metazoa.out", 'w')
##fe = open(args.outprefix+".eukaryota.out", 'w')
##fo = open(args.outprefix+".other.out", 'w')
##
##for seqname in proximity2insecta.keys():
##    ## keys will be seq names
##    total = sum(proximity2insecta[seqname].values())
##    xi = proximity2insecta[seqname]["Insecta"]
##    xm = proximity2insecta[seqname]["Metazoa"]
##    xe = proximity2insecta[seqname]["Eukaryota"]
##    xa = proximity2insecta[seqname]["Archaea"] #will enter this name and 0 into dict if not present
##    xb = proximity2insecta[seqname]["Bacteria"] #will enter this name and 0 into dict if not present
##    proximity2insecta[seqname]["Other"] = total - xi-xm-xe-xa-xb ## Other is a redundant count -- the other thing(s) will be counted by w/e is in that field in above loop
##    outline = ("\t").join([seqname] + [str(k)+":"+str(v) for k,v in proximity2insecta[seqname].iteritems()] + [str(k)+":"+str(v) for k,v in species_counts[seqname].iteritems()])
##    if xi >= args.xi:
##        fi.write(outline+"\n")
##    elif xm >= args.xm:
##        fm.write(outline+"\n")
##    elif xe >= args.xe:
##        fe.write(outline+"\n")
##    else:
##        fo.write(outline+"\n")
##
##            
##fi.close()
##fm.close()
##fe.close()
##fo.close()

