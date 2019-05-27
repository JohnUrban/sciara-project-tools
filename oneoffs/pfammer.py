#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Given interproscan.sh output and either FASTA or GFF,
        add the PFAM description to the FASTA/GFF.

    This utility is only intended to help in post-Maker annotation processing.

    Differences between updates to FASTAs and GFFs:
        - It is assumed one already used Maker's utility ipr_update_gff on the GFF
            - Thus GO IDs were added already
            - Therefore, GO IDs are not added to GFF by this script (but they are to FASTA)
        - Note:
            - This only adds information to the lines with matching IDs
            - Not to children of those IDs, etc
    
    interproscan.sh output interpretation (https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats):
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3. Sequence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)




    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-p', '--pfam', 
                   type= str, 
                   help='''Path to interproscan.sh TSV output.
                        ''')
intype = parser.add_mutually_exclusive_group(required=True)

intype.add_argument('-f', '--fasta', type=str, default=False,
                    help='''Path to single fasta.''')
intype.add_argument('-g', '--gff', type=str, default=False,
                    help='''Path to single gff.''')


args = parser.parse_args()





        
class InterProEntry(object):
    def __init__(self, entry):
        self.entry = entry
    def name(self):
        return self.entry[0].strip()
    def loc(self):
        return self.entry[6].strip() + '-' + self.entry[7].strip()
    def sig(self):
        return self.entry[5].strip()
    def ipa(self):
        return self.entry[12].strip()
    def go(self):
        return ','.join([term.split(':')[1] for term in self.entry[13].strip().split('|')])
    def pa(self):
        return self.entry[14].strip()


class InterPro(object):
    # Will hold dictionary of unique names
    # Each unique name can be associated with multiple entries
    # This will  have functions to extract summary info from each unique name
    def __init__(self):
        self.all = defaultdict(list)
    def add(self, entry):
        entry = entry.split('\t')
        self.all[entry[0]].append( InterProEntry(entry) )
    def get_sig(self, name):
        ans = set([])
        for entry in self.all[name]:
            try:
                ans.add(entry.sig())
            except IndexError:
                pass
        return list(ans)
    def get_ipa(self, name):
        ans = set([])
        for entry in self.all[name]:
            try:
                ans.add(entry.ipa())
            except IndexError:
                pass
        return list(ans)
    def get_go(self, name):
        ans = set([])
        for entry in self.all[name]:
            try:
                ans.add(entry.go())
            except IndexError:
                pass
        return list(ans)
    def get_pa(self, name):
        ans = set([])
        for entry in self.all[name]:
            try:
                ans.add(entry.pa())
            except IndexError:
                pass
        return list(ans)
    def contains(self, name):
        if name in self.all.keys():
            return True
        return False

            
def pfamreader(f):
    d = InterPro()
    with open(f) as pfam:
        for line in pfam:
            d.add(line)
    return d


def update_fa_desc(fa, sfx, key="Note:"):
    if sfx:
        fa.description += '\t'+key+';'.join(sfx)
    return fa

def update_gff_desc(desc, sfx, key="Note:"):
    if sfx:
        desc.append( key+','.join(sfx) )
    return desc

def get_id(desc):
    for e in desc:
        if e.startswith('ID='):
            return e.split('=')[1]
        return None

def process_fasta(fasta, pfam):
    for fa in SeqIO.parse(fasta, 'fasta'):
        if pfam.contains(fa.id):
            fa = update_fa_desc(fa, sfx=pfam.get_sig(fa.id), key="PfamSignature:")
            fa = update_fa_desc(fa, sfx=pfam.get_ipa(fa.id), key="InterProAnnotation:")
            fa = update_fa_desc(fa, sfx=pfam.get_go(fa.id), key="GO:")        
        print ">"+str(fa.description)
        print str(fa.seq)

def process_gff(gff, pfam):
    with open(gff) as f:
        for line in f:
            if line:
                if line.startswith('#'):
                    print line
                else:
                    line = line.strip().split('\t')
                    if len(line)>=9: # i.e. has descriptions column in col9
                        desc = line[8].rstrip(';').split(';')
                        ID = get_id(desc)
                        desc = update_gff_desc(desc, sfx=pfam.get_sig(ID), key="PfamSignature=")
                        desc = update_gff_desc(desc, sfx=pfam.get_ipa(ID), key="InterProAnnotation=")
                    #line = line[:8] + [desc] + line[9:]
                    line[8] = ';'.join(desc)
                    print '\t'.join(line)
                    
            
## Read Pfam
pfam = pfamreader(args.pfam)


try:
    if args.fasta:
        process_fasta(args.fasta, pfam)
        
    elif args.gff:
        process_gff(args.gff, pfam)
except IOError:
    # Broken Pipe
    pass

