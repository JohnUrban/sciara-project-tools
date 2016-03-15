#!/usr/bin/env python

import sys
import argparse
import gzip
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles

parser = argparse.ArgumentParser(description="""

DESCRIPTION
  Takes in tab-separated value (tsv) file with at least 2 fields:
  (1) name of sequence being characterized
  (2) NCBI Tax ID
  Either can be in any column/field of a multi-column tsv file.

  Outputs tsv file with 4 columns
  (1) Name
  (2) Closest classifier to Insecta (Insecta, Metazoa, Eukaryota, Bacteria or Archaea, other)
  (3) species level classification

  This output can be the input for cnv-taxonomy-summarizer.py
    """, formatter_class= argparse.RawTextHelpFormatter)
parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--inputfile', '-i',
                   type=str, default=False,
                   help='''Provide tab-separated file that has taxID in an arbitrary column specified by -k. File can be gzipped''')
parser.add_argument('--taxidcolumn', '-tc',
                   type=int, default=2, required=True,
                   help='''Provide 1-based column number that has the taxIDs in the input file "-i". Default is "-t 2"''')
parser.add_argument('--namecolumn', '-nc',
                   type=int, default=1, required=True,
                   help='''Provide 1-based column number that has the name of sequence being characterized in the input file "-i". Default is "-t 1"''')
parser_input.add_argument('--cmdline', '-c',
                   type=str, default=False,
                   help='''Enter single or comma-separated list of taxIDs at commandline. Cannot be used with "-i".''')
parser.add_argument('--nodes', '-no',
                   type=str, default=False, required=True,
                   help='''Provide path to nodes.dmp or nodes.dmp.gz (from taxdump.tar.gz)''')
parser.add_argument('--names', '-na',
                   type=str, default=False, required=True,
                   help='''Provide path to names.dmp or names.dmp.gz (from taxdump.tar.gz)''')
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

def getTaxTree():
    tree = NcbiTaxonomyFromFiles(open('nodes.dmp'), open('names.dmp'))
    return tree

def getTreeRoot(tree):
    return tree.Root

def taxIDfrom(clade, tree):
    ''' clade e.g. 'Eukaryota' '''
    root = getTreeRoot(tree)
    clade = root.getNodeMatchingName(clade)
    for node in clade.getRankedDescendants('phylum'):
        print node.Name, node.TaxonId

def fullLineageFromTaxID(tree, taxID):
    '''returns lineage information in ranks order'''
    node = tree.ById[taxID]
    lineage = []
    ranks = []
    curr = node
    while curr.Parent is not None:
        lineage.append(curr.Name)
        ranks.append(curr.Rank)
        curr = curr.Parent
    return lineage[-1::-1], ranks[-1::-1]


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
sys.stderr.write("Constructing Taxonomy Tree....\n")
tree = NcbiTaxonomyFromFiles(open(args.nodes), open(args.names))
sys.stderr.write("Done constructing Taxonomy Tree....\n")
if args.cmdline:
    taxIDs = [int(e) for e in args.cmdline.split(",")]
    for taxID in taxIDs:
        taxonomy, ranks = fullLineageFromTaxID(tree, taxID)
        print taxonomy
elif args.inputfile:
    fh = open_file_connection(args.inputfile)
    args.taxidcolumn = args.taxidcolumn - 1
    args.namecolumn = args.namecolumn - 1
    for line in fh:
        line = line.strip().split()
        try:
            taxID = int(line[args.taxidcolumn])
        except ValueError:
            taxID = int(line[args.taxidcolumn].split(";")[0])
        name = line[args.namecolumn]
        try:
            taxonomy, ranks = fullLineageFromTaxID(tree, taxID)
        except KeyError:
            ## if anything raises an error, move on to next
            sys.stderr.write(name +" had a taxID raised an KeyError: " + str(taxID)+"\n")
            continue
        species = ("_").join(taxonomy[-1].split())
        if "Insecta" in taxonomy:
            print ("\t").join([name, "Insecta", species])
        elif "Metazoa" in taxonomy:
            print ("\t").join([name, "Metazoa", species])
        elif "Eukaryota" in taxonomy:
            print ("\t").join([name, "Eukaryota", species])
        elif "Bacteria" in taxonomy:
            print ("\t").join([name, "Bacteria", species])
        elif "Archaea" in taxonomy:
            print ("\t").join([name, "Archaea", species])
        else:
            print ("\t").join([name] + taxonomy[:2])
            


