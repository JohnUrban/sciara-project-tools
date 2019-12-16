#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Given tables and either FASTA or GFF,
        update the FASTA/GFF.

    This utility is only intended to help in post-Maker annotation processing.

    Also see pfammer.py



    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-m', '--map', 
                   type= str, 
                   help='''Path to interproscan.sh TSV output.
                        ''')
parser.add_argument('-m1', '--mapcol1', 
                   type=int, default=0,
                   help='''0-based column for key (that will be updated in the FA or GFF).
                        ''')
parser.add_argument('-m2', '--mapcol2', 
                   type= str, 
                   help='''0-based column for paired value (description).
                        ''')

parser.add_argument('-k', '--key', 
                   type=str, required=True,
                   help='''What is added to GFF: key=stuff.
                           What is added to Fasta: key:stuff.
                        ''')

intype = parser.add_mutually_exclusive_group(required=True)

intype.add_argument('-f', '--fasta', type=str, default=False,
                    help='''Path to single fasta.''')
intype.add_argument('-g', '--gff', type=str, default=False,
                    help='''Path to single gff.''')


args = parser.parse_args()





class Map(object):
    # Will hold dictionary of unique names
    # Each unique name can be associated with multiple entries
    # This will  have functions to extract summary info from each unique name
    def __init__(self, c1=0, c2=1):
        self.c1 = c1
        self.c2 = c2
        self.all = defaultdict(list)
    def add(self, entry):
        entry = entry.strip().split('\t')
        self.all[entry[self.c1]].append( entry[self.c2] )
    def contains(self, name):
        if name in self.all.keys():
            return True
        return False
    def get(self, key):
        try:
            return self.all[key]
        except:
            return False

            
def mapreader(f):
    d = Map()
    with open(f) as mapfile:
        for line in mapfile:
            d.add(line)
    return d


def update_fa_desc(fa, update, key="Note:"):
    if update:
        fa.description += '\t'+key+';'.join(update)
    return fa

def update_gff_desc(desc, update, key="Note:"):
    if update:
        desc.append( key+','.join(update) )
    return desc

def get_id(desc):
    for e in desc:
        if e.startswith('ID='):
            return e.split('=')[1]
        return None

def process_fasta(fasta, observed, key="Note"):
    for fa in SeqIO.parse(fasta, 'fasta'):
        if observed.contains(fa.id):
            fa = update_fa_desc(fa, update=observed.get(fa.id), key=key+":")   
        print ">"+str(fa.description)
        print str(fa.seq)

def process_gff(gff, observed, key="Note"):
    with open(gff) as f:
        for line in f:
            if line:
                if line.startswith('#'):
                    print line.rstrip()
                else:
                    line = line.strip().split('\t')
                    if len(line)>=9: # i.e. has descriptions column in col9
                        desc = line[8].rstrip(';').split(';')
                        ID = get_id(desc)
                        desc = update_gff_desc(desc, update=observed.get(ID), key=key+":")
                    #line = line[:8] + [desc] + line[9:]
                    line[8] = ';'.join(desc)
                    print '\t'.join(line)
                    
            
## Read Map
observed = mapreader(args.map)


try:
    if args.fasta:
        process_fasta(args.fasta, observed, key=args.key)
        
    elif args.gff:
        process_gff(args.gff, observed, key=args.key)
except IOError:
    # Broken Pipe
    pass

