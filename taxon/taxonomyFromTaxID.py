#!/usr/bin/env python

import sys
import argparse
import gzip
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles
import cPickle as pickle

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


  TODO:
  Also include bitscore information, and use that in summarizing.
  Allow option to not only add counts (and bitscores) to closest taxonomy level, but to everything above it too.
  This would allow different minimum level counts. e.g., if order Diptera is found (+1) then
  class Insecta, phylum arthropoda, kingdom metazoa, and super kingdom Eukaryota would all get +1 as well.
  This way, it can be more sensitive to the 'closest hit' while allowing more stringency, since if
  Diptera is given +1 under the current approach, all parental levels are not changed - leaving the following possibility
  among others.
  Assume you have 7 hits for a sequence, with the following closest levels:
  1 Diptera, 1 Insecta, 1 Arthropoda, 1 Metazoa, 1 Eukaryota, 2 Bacteria.
  And say you set the scoring scheme to requiring a minimum of 2 hits for all levels, including other super kingdoms.
  Then the sequence would be summarized as bacterial.
  Under the parental-inheritance of hits scheme, the counts would look like:
  1 Diptera, 2 Insecta, 3 Arthropoda, 4 Metazoa, 5 Eukaryota, 2 Bacteria.
  With the same min levels scoring scheme, the sequence would be labeled as Insecta
  (or Arthropoda if one wanted a tie-breaker).
  Bitscores (or expected valies) will help as well.
  The same hits above might have the following bit scores:
  Diptera 30, Insect 34, Arthropoda 56, Metazoa 41, Eukaryota 59, Bacteria 88 and 112
  Using the highest bitscore, Bacteria would win.
  Using sum of bit scores with parental inheritance:
  Diptera 30, Insecta 64, Arthropoda 120, Metazoa 161, Eukaryota 220, Bacteria 200
  Eukaroyota would win --- though I do not know whether summing or single best makes most sense.

  Could also include length of hit. Or coordinates along seq of hits to see if the seq is a missassembly...
  or for further experiments to see if it is HGT.

  Allow TaxIDs to be used for levels.
  Allow species TaxID to be submitted alone with option to define all other levels automatically using the TaxTree.
  
  John Urban (2015, 2016)
    """, formatter_class= argparse.RawTextHelpFormatter)
parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--inputfile', '-i',
                   type=str, default=False,
                   help='''Provide tab-separated file that has taxID in an arbitrary column specified by -k. File can be gzipped''')
parser_input.add_argument('--pickleonly', '-p', type=str, default=False,
                    help=''' Make taxonomy tree, pickle it, save in given file name(to also use later with --taxtree), then exit. Note: picking and loading does not make it faster.''')
parser.add_argument('--pickle', type=str, default=False,
                    help=''' Make taxonomy tree, pickle it, save in given file name (to also use later with --taxtree). But don't exit. Continue processing inputfile.''')
parser.add_argument('--taxidcolumn', '-tc',
                   type=int, default=2, required=False,
                   help='''Provide 1-based column number that has the taxIDs in the input file "-i". Default is "-t 2"''')
parser.add_argument('--namecolumn', '-nc',
                   type=int, default=1, required=False,
                   help='''Provide 1-based column number that has the name of sequence being characterized in the input file "-i". Default is "-t 1"''')
parser_input.add_argument('--cmdline', '-c',
                   type=str, default=False,
                   help='''Enter single or comma-separated list of taxIDs at commandline. Cannot be used with "-i".''')
parser.add_argument('--nodes', '-no',
                   type=str, default=False, required=False,
                   help='''Provide path to nodes.dmp or nodes.dmp.gz (from taxdump.tar.gz)''')
parser.add_argument('--names', '-na',
                   type=str, default=False, required=False,
                   help='''Provide path to names.dmp or names.dmp.gz (from taxdump.tar.gz)''')
parser.add_argument('--taxtree', type=str, default=False, required=False,
                    help='''Instead of computing taxonomy tree from names/nodes, load in a pre-computed, pickled tree.''')

parser.add_argument('--superkingdom', type=str, default="Eukaryota",
                    help='''Provide super kingdom of your target species. Default: Eukaryota. Example of another: Bacteria, Archaea.''')

parser.add_argument('--kingdom', type=str, default="Metazoa",
                    help='''Provide kingdom of your target species. Default: Metazoa.''')


parser.add_argument('--phylum', type=str, default="Arthropoda",
                    help='''Provide phylum of your target species. Default: Eukaryota.
                    Example of others: Apicomplexa, Arthropoda, Ascomycota, Basidiomycota, Blastocladiomycota, Chordata,
                    Cnidaria, Echinodermata, Eukaryota-undef, Glomeromycota, Mollusca, Nematoda, Platyhelminthes,
                    Porifera, Proteobacteria, Streptophyta, Xenacoelomorpha''')

parser.add_argument('--Class', type=str, default="Insecta",
                    help='''Provide class of your target species. Default: Insecta.''')


parser.add_argument('--order', type=str, default="Diptera",
                    help='''Provide order of your target species. Default: Diptera.
                    Example of others:
                    Astigmata, Beloniformes, Chiroptera, Coleoptera, Dasyuromorphia, Decapoda, Echinoida,
                    Enterogona, Fabales, Glomerales,Hemiptera, Hymenoptera, Insectivora, Lepidoptera,
                    Mesostigmata, Octopoda, Phthiraptera, Plecoptera, Primates, Rickettsiales, Rodentia,
                    Saccharomycetales, Thysanoptera, Xiphosura''')


parser.add_argument('--family', type=str, default="Sciaridae",
                    help='''Provide order of your target species. Default: Sciaridae.
                    Example of others:
                    Aphididae, Apidae, Apterygidae, Bombycidae, Bombyliidae, Bovidae, Chironomidae,
                    Cichlidae, Diopsidae, Drosophilidae, Eurypygidae, Formicidae, Glomeraceae, Hominidae,
                    Muscidae, Nymphalidae, Octopodidae, Pachyneuridae, Rickettsiaceae, Saccharomycetaceae,
                    Trichogrammatidae, Vespidae, Xylophagidae''')

parser.add_argument('--genus', type=str, default="Bradysia",
                    help='''Provide order of your target species. Default: Bradysia.
                    Example of others:
                    Aedes, Anopheles, Aphidius, Apis, Bombus, Bombyx, Bradysia, Chironomus, Drosophila, Homo,
                    Mus, Musca, Nasonia, Nicrophorus, Nilaparvata, Rhynchosciara, Saccharomyces, Sciara''')

parser.add_argument('--species', type=str, default="Bradysia coprophila",
                    help='''Provide order of your target species. Default: Bradysia coprophila.
                    Example of others:
                    Aedes aegypti, Anopheles gambiae, Apis mellifera, Arabidopsis thaliana, Bombyx mori, Drosophila melanogaster, Drosophila miranda...''')

parser.add_argument('--levels', type=str, default='species,genus,family,order,class,phylum,kingdom,superkingdom,othersuperkingdoms',
                    help='''Use only the given levels for classification. Give comma-separated list.
                    Default (all):
                    species,genus,family,order,class,phylum,kingdom,superkingdom,othersuperkingdoms''')

args = parser.parse_args()

assert args.taxtree or (args.nodes and args.names)
assert (args.taxtree and not (args.nodes and args.names)) or (not args.taxtree and (args.nodes and args.names))
if args.inputfile:
    assert args.taxidcolumn and args.namecolumn
    ## since have args.pickleonly option, these are not "required"; but it is safer to require them when inputfile used.

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

if args.names and args.nodes:
    sys.stderr.write("Constructing Taxonomy Tree....\n")
    tree = NcbiTaxonomyFromFiles(open(args.nodes), open(args.names))
    sys.stderr.write("Done constructing Taxonomy Tree....\n")
    if args.pickle or args.pickleonly:
        sys.stderr.write("Pickling Taxonomy Tree....\n")
        if args.pickle:
            out_fname = args.pickle
        else:
            out_fname = args.pickleonly
        with open(out_fname, 'wb') as out:
            pickle.dump(tree, out, -1)
        sys.stderr.write("Done pickling Taxonomy Tree....\n")
        if args.pickleonly:
            quit()
elif args.taxtree:
    sys.stderr.write("Loading Taxonomy Tree....\n")
    with open(args.taxtree, 'rb') as pkl_file:
        tree = pickle.load(pkl_file)
    sys.stderr.write("Done loading Taxonomy Tree....\n")


all_levels = 'species,genus,family,order,class,phylum,kingdom,superkingdom'.split(',')
all_level_args = [args.species, args.genus, args.family, args.order, args.Class, args.phylum, args.kingdom, args.superkingdom]
levels = []
for i in range(len(all_levels)):
    if all_levels[i] in args.levels:
        levels.append(all_level_args[i])
if "othersuperkingdoms" in args.levels:
    other_super_kingdoms = ["Eukaryota", "Bacteria", "Archaea"]
    other_super_kingdoms.remove(args.superkingdom)
    levels += other_super_kingdoms


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
##            print taxonomy, ranks
        except KeyError:
            ## if anything raises an error, move on to next
            sys.stderr.write(name +" had a taxID raised an KeyError: " + str(taxID)+"\n")
            continue
        try:
            species = ("_").join(taxonomy[-1].split())
        except:
            species = "NA"
        level_found = False
        for level in levels:
            if level in taxonomy:
                print ("\t").join([name, level, species])
                level_found = True
                break
        if not level_found:
            try:
                sk = taxonomy[1]
            except:
                sk = "NA"
            print ("\t").join([name] + [sk, species])
                
                
        
##        if "Insecta" in taxonomy:
##            print ("\t").join([name, "Insecta", species])
##        elif "Metazoa" in taxonomy:
##            print ("\t").join([name, "Metazoa", species])
##        elif args.superkingdom in taxonomy:
##            print ("\t").join([name, "Eukaryota", species])
##        elif "Bacteria" in taxonomy:
##            print ("\t").join([name, "Bacteria", species])
##        elif "Archaea" in taxonomy:
##            print ("\t").join([name, "Archaea", species])
##        else:
##            print ("\t").join([name] + taxonomy[:2])
            


