#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict
import numpy as np

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


Maker QI columns
1       =       Length of the 5 UTR
2       =       Fraction of splice sites confirmed by an EST alignment
3       =       Fraction of exons that overlap an EST alignment
4       =       Fraction of exons that overlap EST or Protein alignments
5       =       Fraction of splice sites confirmed by a SNAP prediction
6       =       Fraction of exons that overlap a SNAP prediction
7       =       Number of exons in the mRNA
8       =       Length of the 3 UTR
9       =       Length of the protein sequence produced by the mRNA

    """, formatter_class= argparse.RawTextHelpFormatter)



intype = parser.add_mutually_exclusive_group(required=True)

intype.add_argument('-m', '--maker', type=str, default=False,
                    help='''Path to Maker GFF.''')
intype.add_argument('-s', '--stringtie', type=str, default=False,
                    help='''Path to StringTie GFF.''')

parser.add_argument('-P', '--prefix',
                   type=str, default=False,
                   help='''If prefix is provided, single-column text files will be output with the given prefix.
                        These will include (where applicable):
                            Number of exons per gene
                            Gene lengths
                            RNA lengths
                            Exon lengths
                            Intron lengths
                            5' UTR lengths
                            3' UTR lengths
                        ''')
parser.add_argument('-F', '--full',
                   action='store_true', default=False,
                   help='''--maker is full by default, meaning it outputs info on CDS, 5pUTR, 3pUTR, ontology, pfam, blastp, etc...
                        --stringtie assumes it is not full.
                        However, there is nothing stopping one from updating StringTie GFFs to have this information.
                        For the moment, the StringTie option would not interpret much of it... TODO.
                        The recommendation is if you plan to update a StringTie GFF with that information,
                            then fully convert it into a Maker-like GFF and use the --maker option.
                        ''')

args = parser.parse_args()



class GFF_Entry(object):
    def __init__(self, entry, attributes):
        ''' entry       = GFF line split on tabs
            attributes  = dict output from GFF._parse_desc()

            The following attributes are anticipated though may not be present:
                Alias
                Dbxref
                ID
                Name
                Note
                Ontology_term
        '''
        self.entry = entry
        self.attributes_dict = attributes
        self.has_attribute = {}
    def seqid(self):
        return self.entry[0]
    def source(self):
        return self.entry[1]
    def type(self):
        return self.entry[2]
    def start(self):
        return int(self.entry[3])
    def end(self):
        return int(self.entry[4])
    def score(self):
        return self.entry[5]
    def strand(self):
        return self.entry[6]
    def phase(self):
        return self.entry[7]
    def attributes(self):
        return self.entry[8]
    def get_attribute(self, key):
        if self.in_attributes(key):
            return self.attributes_dict[key]
    def in_attributes(self, key):
        try:
            return self.has_attribute[key]
        except:
            self.has_attribute[key] = key in self.attributes_dict.keys()
            return self.has_attribute[key]

    def attribute_is(self, key, value):
        return self.get_attribute(key) == value
    
    def ID(self):
        '''All expected to have ID'''
        return self.get_attribute('ID')
    def length(self):
        return self.end() - self.start() + 1 #self.end() - self.start + 1

    


class Gene(GFF_Entry):
    def __init__(self, entry, attributes):
        ''' entry       = GFF line split on tabs
            attributes  = dict output from GFF._parse_desc()

            The following attributes are anticipated though may not be present:
                Alias
                Dbxref
                ID
                Name
                Note
                Ontology_term
        '''
        GFF_Entry.__init__(self, entry, attributes)
        self.mRNA = {}
        self.exons = {}
        self.CDS = {}
        self.fiveUTR = {}
        self.threeUTR = {}
        self.introns = None
        self.intron_lengths = None

    def add_mRNA(self, RNA):
        ''' RNA     =   mRNA() object'''
        self.mRNA[RNA.ID()] = RNA
    def add_exon(self, exon):
        ''' RNA     =   mRNA() object'''
        self.exons[exon.ID()] = exon
    def add_CDS(self, cds):
        ''' RNA     =   mRNA() object'''
        self.CDS[cds.ID()] = cds
    def add_5pUTR(self, utr):
        ''' RNA     =   mRNA() object'''
        self.fiveUTR[utr.ID()] = utr
    def add_3pUTR(self, utr):
        ''' RNA     =   mRNA() object'''
        self.threeUTR[utr.ID()] = utr

    def num_mRNA(self):
        return len(self.mRNA.keys())
    def num_exons(self):
        return len(self.exons.keys())
    def num_cds(self):
        return len(self.CDS.keys())
    def num_5pUTR(self):
        return len(self.fiveUTR.keys())
    def num_3pUTR(self):
        return len(self.threeUTR.keys())
    def define_merged_exons(self):
        exons = []
        self.introns = []
        self.intron_lengths = []
        for exon in self.exons.keys():
            exons.append( [self.exons[exon].start(), self.exons[exon].end()] )
        exons = sorted(exons)
        merged_exons = []
        curr_exon = exons[0]
        for i in range(1, len(exons)):
            old_exon = curr_exon
            curr_exon = exons[i]
            if curr_exon[0] < old_exon[1] and curr_exon[0] >= old_exon[0]:
                if curr_exon[1] > old_exon[1]:
                    curr_exon = [old_exon[0], curr_exon[1]]
                else:
                    curr_exon = old_exon
            else:
                merged_exons.append( old_exon )
        merged_exons.append( curr_exon )
        return sorted(merged_exons)
                
    def define_introns(self):
##        exons = []
##        for exon in self.exons.keys():
##            exons.append( (self.exons[exon].start(), self.exons[exon].end()) )
##        exons = sorted(exons)
        exons = self.define_merged_exons()
        for i in range(1, len(exons)):
            start = exons[i-1][1]+1
            end = exons[i][0]-1
            self.introns.append( (start, end) )
            self.intron_lengths.append( end - start + 1 )

            ##DEBUG
##            if end - start + 1 <= 0:
##                print "ERROR:"
##                print exons
##                print mexons
##                print start, end
##                print end - start + 1
##                print 
            
    def num_introns(self):
        if self.introns is None:
            self.define_introns()
        return len(self.introns)

    def get_intron_lengths(self):
        if self.introns is None:
            self.define_introns()
        return self.intron_lengths
            
    def print_exon_intron_structure(self):
        ## TODO: exons in gene can overlap etc, so not totally correct here
        ##  Can merge all -- also not totally correct since it is excluding the "0 bp introns"
##        exons = []
##        for exon in self.exons.keys():
##            exons.append( (self.exons[exon].start(), self.exons[exon].end()) )
##        exons = sorted(exons)
        print "Gene", self.start(), self.end()
        exons = self.define_merged_exons()
        print "Exon", exons[0][0], exons[0][1]
        for i in range(1, len(exons)):
            print "Intron", exons[i-1][1]+1, exons[i][0]-1
            print "Exon", exons[i][0], exons[i][1]

    def in_mRNA_attrbibutes(self, key):
        mRNAhas = False
        for ID in self.mRNA.keys():
            if self.mRNA[ID].in_attributes(key):
                mRNAhas = True
        return mRNAhas

    def in_mRNA_attrbibutes_not_counting_given_value(self, key, value):
        mRNAhas = False
        for ID in self.mRNA.keys():
            if self.mRNA[ID].in_attributes(key) and not self.mRNA[ID].attribute_is(key, value):
                mRNAhas = True
        return mRNAhas

    def has(self, key):
        return self.in_attributes(key) or self.in_mRNA_attrbibutes(key)

    def get_mRNA_attribute(self, key):
        mRNAhas = []
        for ID in self.mRNA.keys():
            mRNAhas.append( self.mRNA[ID].get_attribute(key) )
        return mRNAhas
    
    def has_but_not_counting_given_value(self, key, value):
        genehas = self.has(key) and not self.attribute_is(key, value)
        mRNAhas = self.in_mRNA_attrbibutes_not_counting_given_value(key, value)
        #print genehas, mRNAhas, key, value, self.attribute_is(key, value), self.get_attribute(key), self.get_mRNA_attribute(key)
        return genehas or mRNAhas
        
    def has_ontology_term(self):
        return self.has('Ontology_term')

    def has_blastp(self):
        return self.has_but_not_counting_given_value('Note', 'Protein of unknown function')

    def has_pfam(self):
        return self.has('Dbxref') or self.has('PfamSignature') or self.has('InterProAnnotation')
          
    def has_ontology_blastp_and_pfam(self):
        return self.has_ontology_term() and self.has_blastp() and self.has_pfam()

    def has_ontology_blastp_or_pfam(self):
        return self.has_ontology_term() or self.has_blastp() or self.has_pfam()

    def update_gene_coords(self):
        ''' This is mainly for the StringTieGFF class.
            StringTie GFFs only have transcripts and exons.
            However, transcript attributes have a geneID, which is analogous to the mRNA's Parent gene for Maker.
            Thus, the Gene object is initialized with the first transcript and updated with each transcript added'''
        start = self.start()
        end = self.end()
        for ID in self.mRNA.keys():
            if self.mRNA[ID].start() < start:
                start = self.mRNA[ID].start()
            if self.mRNA[ID].end() > end:
                start = self.mRNA[ID].end()
        self.entry[3] = str(start)
        self.entry[4] = str(end)

class GenePart(GFF_Entry):
    def __init__(self, entry, attributes):
        ''' entry       = GFF line split on tabs
            attributes  = dict output from GFF._parse_desc()

            Should be used for GFF lines types where the following attributes
            can be anticipated though may not be present:
                ID
                Parent

            For Maker, the types are:
                exon
                CDS
                five_prime_UTR
                three_prime_UTR

        '''
        GFF_Entry.__init__(self, entry, attributes)
    def Parent(self):
        return self.get_attribute('Parent')
    def Parents(self, delim=','):
        return self.get_attribute('Parent').split(delim)


            
class mRNA(GenePart):
    def __init__(self, entry, attributes):
        '''
            entry       = GFF line split on tabs
            attributes  = dict output from GFF._parse_desc()

            The following attributes are anticipated though may not be present:
                Alias
                Dbxref
                ID
                Name
                Note
                Ontology_term
                Parent
                _AED
                _QI
                _eAED
                PfamSignature
                InterProAnnotation


            The QI columns from Maker are:
                0       =       Length of the 5 UTR
                1       =       Fraction of splice sites confirmed by an EST alignment
                2       =       Fraction of exons that overlap an EST alignment
                3       =       Fraction of exons that overlap EST or Protein alignments
                4       =       Fraction of splice sites confirmed by a SNAP prediction
                5       =       Fraction of exons that overlap a SNAP prediction
                6       =       Number of exons in the mRNA
                7       =       Length of the 3 UTR
                8       =       Length of the protein sequence produced by the mRNA
        '''
        GFF_Entry.__init__(self, entry, attributes)
        self.QIlist = None
        
    def Alias(self):
        return self.get_attribute('Alias')
    def Dbxref(self):
        return self.get_attribute('Dbxref')
    def Name(self):
        return self.get_attribute('Name')
    def Note(self):
        return self.get_attribute('Note')
    def Ontology_term(self):
        return self.get_attribute('Ontology_term')
    def AED(self):
        return float(self.get_attribute('_AED'))
    def QI(self):
        return self.get_attribute('_QI')
    def eAED(self):
        return float(self.get_attribute('_eAED'))
    def PfamSignature(self):
        return self.get_attribute('PfamSignature')
    def InterProAnnotation(self):
        return self.get_attribute('InterProAnnotation')
        
    def get_from_QI(self, idx):
        if self.QIlist is None:
            self.QIlist = self.QI().split('|')
        return self.QIlist[idx]
        
    def get_5pUTR_length_from_QI(self):
        return int(self.get_from_QI(idx=0))
    def get_3pUTR_length_from_QI(self):
        return int(self.get_from_QI(idx=7))
    def get_protein_length_from_QI(self):
        return int(self.get_from_QI(idx=8))
    
    def get_fraction_splice_sites_confirmed_by_EST_from_QI(self):
        return float(self.get_from_QI(idx=1))
    def get_fraction_splice_sites_confirmed_by_SNAP_from_QI(self):
        return float(self.get_from_QI(idx=4))
    
    def get_fraction_exons_that_overlap_EST_from_QI(self):
        return float(self.get_from_QI(idx=2))
    def get_fraction_exons_that_overlap_EST_or_protaln_from_QI(self):
        return float(self.get_from_QI(idx=3))
    def get_fraction_exons_that_overlap_SNAP_from_QI(self):
        return float(self.get_from_QI(idx=5))

    def get_number_of_exons_from_QI(self):
        return int(self.get_from_QI(idx=6))


    def has_ontology_term(self):
        return self.in_attributes('Ontology_term')
    def has_blastp(self):
        return self.in_attributes('Note')

    def has_blastp(self):
        return self.in_attributes('Note') and not self.attribute_is('Note', 'Protein of unknown function')
    
    def has_pfam(self):
        return self.in_attributes('Dbxref') or self.in_attributes('PfamSignature') or self.in_attributes('InterProAnnotation')

    def has_ontology_blastp_and_pfam(self):
        return self.has_ontology_term() and self.has_blastp() and self.has_pfam()

    def has_ontology_blastp_or_pfam(self):
        return self.has_ontology_term() or self.has_blastp() or self.has_pfam()





class GFF(object):
    def __init__(self, allowable_sources=['maker'], allowable_types=['gene', 'mRNA', 'five_prime_UTR', 'exon', 'CDS', 'three_prime_UTR']):
        self.genes = {}
        self.header = ''
        self.headerpassed = False
        self.allowable_sources = allowable_sources
        self.allowable_types = allowable_types
        self.num_genes = 0
        self.num_mRNA = 0
        self.num_exons = 0
        self.gene_lengths = []
        self.mRNA_lengths = []
        self.exon_lengths = []
        self.num_5pUTR = 0
        self.five_prime_UTR_lengths = []
        self.num_3pUTR = 0
        self.three_prime_UTR_lengths = []
        self.intron_lengths = []
        self.num_cds = 0
        self.cds_lengths = []
        self.RNA_parents = defaultdict(set)
        self.exons_per_gene = None
        self.mRNA_per_gene = None
        self.num_introns = None
        self.intron_lengths = None

    def add(self, line):
        ''' line = unsplit GFF line '''
        entry = line.strip().split('\t')
        if line:
            if line.startswith('#') and not self.headerpassed:
                self.header += line
            elif line.startswith('#') and self.headerpassed:
                pass #continue
            elif len(entry)>=9 and entry[1] in self.allowable_sources and entry[2] in self.allowable_types:
                # i.e. has descriptions column in col9
                desc = entry[8].rstrip(';').split(';')
                
                d = self._parse_desc(desc)
                if entry[2] == 'gene':
                    gene = Gene(entry, attributes=d)
                    self.genes[d['ID']] = gene
                    self.num_genes += 1
                    self.gene_lengths.append ( gene.length() )
                elif entry[2] == 'mRNA':
                    RNA = mRNA(entry, attributes=d)
                    self.num_mRNA += 1
                    self.mRNA_lengths.append( RNA.length() )
                    for gene in RNA.Parents():
                        self.genes[gene].add_mRNA(RNA)
                        self.RNA_parents[RNA.ID()].add( gene )
                elif entry[2] == 'exon':
                    exon = GenePart(entry, attributes=d)
                    self.num_exons += 1
                    self.exon_lengths.append( exon.length() )
                    genes = set([])
                    for RNA in exon.Parents():
                        for gene in self.RNA_parents[RNA]:
                            genes.add(gene)
                    for gene in list(genes): #should only be 1
                        self.genes[gene].add_exon( exon )
                elif entry[2] == 'CDS':
                    cds = GenePart(entry, attributes=d)
                    self.num_cds += 1
                    self.cds_lengths.append( cds.length() )
                    genes = set([])
                    for RNA in cds.Parents():
                        for gene in self.RNA_parents[RNA]:
                            genes.add(gene)

                    for gene in list(genes): #should only be 1
                        self.genes[gene].add_CDS( cds )
                elif entry[2] == 'five_prime_UTR':
                    genepart = GenePart(entry, attributes=d)
                    self.num_5pUTR += 1
                    self.five_prime_UTR_lengths.append( genepart.length() )
                    genes = set([])
                    for RNA in genepart.Parents():
                        for gene in self.RNA_parents[RNA]:
                            genes.add(gene)
                    for gene in list(genes): #should only be 1
                        self.genes[gene].add_5pUTR( genepart )
                elif entry[2] == 'three_prime_UTR':
                    genepart = GenePart(entry, attributes=d)
                    self.num_3pUTR += 1
                    self.three_prime_UTR_lengths.append( genepart.length() )
                    genes = set([])
                    for RNA in genepart.Parents():
                        for gene in self.RNA_parents[RNA]:
                            genes.add(gene)
                    for gene in list(genes): #should only be 1
                        self.genes[gene].add_3pUTR( genepart )
    
    def get_num_genes(self):
        return self.num_genes

    def get_gene_lengths(self):
        return self.gene_lengths()

    def get_num_exons(self):
        return self.num_exons

    def get_exon_lengths(self):
        return self.exon_lengths

    def get_num_cds(self):
        return self.num_cds

    def get_cds_lengths(self):
        return self.cds_lengths

    def get_num_mRNA(self):
        return self.num_mRNA

    def get_num_genes_with_ontology(self):
        return sum([self.genes[e].has_ontology_term() for e in self.genes.keys()])

    def get_num_genes_with_blastp(self):
        return sum([self.genes[e].has_blastp() for e in self.genes.keys()])

    def get_num_genes_with_pfam(self):
        return sum([self.genes[e].has_pfam() for e in self.genes.keys()])

    def get_num_genes_with_ontology_blastp_and_pfam(self):
        return sum([self.genes[e].has_ontology_blastp_and_pfam() for e in self.genes.keys()])

    def get_num_genes_with_ontology_blastp_or_pfam(self):
        return sum([self.genes[e].has_ontology_blastp_or_pfam() for e in self.genes.keys()])


    def get_mRNA_lengths(self):
        return self.mRNA_lengths

    def get_num_5pUTR(self):
        return self.num_5pUTR

    def get_5pUTR_lengths(self):
        return self.five_prime_UTR_lengths

    def get_num_3pUTR(self):
        return self.num_3pUTR

    def define_introns(self):
        self.intron_lengths = []
        self.num_introns = 0
        for gene in self.genes.keys():
            self.num_introns += self.genes[gene].num_introns()
            self.intron_lengths += self.genes[gene].get_intron_lengths()

    def get_num_introns(self):
        if self.num_introns is None:
            self.define_introns()
        return self.num_introns
    

    def get_3pUTR_lengths(self):
        return self.three_prime_UTR_lengths

    def get_exons_per_gene(self):
        if self.exons_per_gene is None:
            self.exons_per_gene = [self.genes[gene].num_exons() for gene in self.genes.keys()]
        return self.exons_per_gene

    def get_mRNA_per_gene(self):
        if self.mRNA_per_gene is None:
            self.mRNA_per_gene = [self.genes[gene].num_mRNA() for gene in self.genes.keys()]
        return self.mRNA_per_gene

    def get_stats(self, x, p=[0,0.05, 0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]):
        try:
            q = np.quantile(x, p)
            return [np.mean(x), np.std(x, ddof=1), np.min(x), np.max(x), 'Q:'] + [e for e in  q]
        except:
            return None
    def get_exon_per_gene_stats(self):
        return self.get_stats(np.array(self.get_exons_per_gene()))

    def get_mRNA_per_gene_stats(self):
        return self.get_stats(np.array(self.get_mRNA_per_gene()))
        

    def get_exon_length_stats(self):
        return self.get_stats(np.array(self.exon_lengths))
        
    def get_cds_length_stats(self):
        return self.get_stats(np.array(self.cds_lengths))

    def get_mRNA_length_stats(self):
        return self.get_stats(np.array(self.mRNA_lengths))

    def get_gene_length_stats(self):
        if not self.gene_lengths:
            self.gene_lengths = [self.genes[gene].length() for gene in self.genes.keys()]
        return self.get_stats(np.array(self.gene_lengths))

    def get_5pUTR_length_stats(self):
        return self.get_stats(np.array(self.five_prime_UTR_lengths))

    def get_3pUTR_length_stats(self):
        return self.get_stats(np.array(self.three_prime_UTR_lengths))

    def get_intron_length_stats(self):
        if self.intron_lengths is None:
            self.define_introns()
        return self.get_stats(np.array(self.intron_lengths))

    def get_gene_lengths(self):
        if not self.gene_lengths:
            self.gene_lengths = [self.genes[gene].length() for gene in self.genes.keys()]
        return self.gene_lengths

    def get_mRNA_lengths(self):
        return self.mRNA_lengths
    
    def get_exon_lengths(self):
        return self.exon_lengths

    def get_5pUTR_lengths(self):
        return self.five_prime_UTR_lengths

    def get_3pUTR_lengths(self):
        return self.three_prime_UTR_lengths

    def get_intron_lengths(self):
        if self.intron_lengths is None:
            self.define_introns()
        return self.intron_lengths
    

    def print_exon_intron_structure(self):
        for gene in self.genes.keys():
            self.genes[gene].print_exon_intron_structure()


    def get_num_5pUTR_per_gene(self):
        return [self.genes[gene].num_5pUTR() for gene in self.genes.keys()]

    def get_num_3pUTR_per_gene(self):
        return [self.genes[gene].num_3pUTR() for gene in self.genes.keys()]

    def get_num_5pUTR_per_gene_stats(self):
        return self.get_stats(np.array(self.get_num_5pUTR_per_gene()))

    def get_num_3pUTR_per_gene_stats(self):
        return self.get_stats(np.array(self.get_num_3pUTR_per_gene()))

    def get_num_genes_with_5pUTR(self):
        return sum([self.genes[gene].num_5pUTR() > 0 for gene in self.genes.keys()])
    def get_num_genes_with_3pUTR(self):
        return sum([self.genes[gene].num_3pUTR() > 0 for gene in self.genes.keys()])


    
    def get_fraction_splice_sites_confirmed_by_EST_from_QI(self):
        return [self.genes[geneID].mRNA[rnaID].get_fraction_splice_sites_confirmed_by_EST_from_QI() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]
        
    def get_fraction_splice_sites_confirmed_by_SNAP_from_QI(self):
        return [self.genes[geneID].mRNA[rnaID].get_fraction_splice_sites_confirmed_by_SNAP_from_QI() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]
        
    def get_fraction_exons_that_overlap_EST_from_QI(self):
        return [self.genes[geneID].mRNA[rnaID].get_fraction_exons_that_overlap_EST_from_QI() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]    

    def get_fraction_exons_that_overlap_EST_or_protaln_from_QI(self):
        return [self.genes[geneID].mRNA[rnaID].get_fraction_exons_that_overlap_EST_or_protaln_from_QI() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]    

    def get_fraction_exons_that_overlap_SNAP_from_QI(self):
        return [self.genes[geneID].mRNA[rnaID].get_fraction_exons_that_overlap_SNAP_from_QI() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]

    def get_AEDs(self):
        return [self.genes[geneID].mRNA[rnaID].AED() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]

    def get_eAEDs(self):
        return [self.genes[geneID].mRNA[rnaID].eAED() for geneID in self.genes.keys() for rnaID in self.genes[geneID].mRNA.keys()]
    

    def get_fraction_splice_sites_confirmed_by_EST_from_QI_stats(self):
        return self.get_stats(np.array( self.get_fraction_splice_sites_confirmed_by_EST_from_QI() ) )
    
        
    def get_fraction_splice_sites_confirmed_by_SNAP_from_QI_stats(self):
        return self.get_stats(np.array( self.get_fraction_splice_sites_confirmed_by_SNAP_from_QI() ) )
    
    def get_fraction_exons_that_overlap_EST_from_QI_stats(self):
        return self.get_stats(np.array( self.get_fraction_exons_that_overlap_EST_from_QI() ) )
                              
    def get_fraction_exons_that_overlap_EST_or_protaln_from_QI_stats(self):
        return self.get_stats(np.array( self.get_fraction_exons_that_overlap_EST_or_protaln_from_QI() ) )
                              
    def get_fraction_exons_that_overlap_SNAP_from_QI_stats(self):
        return self.get_stats(np.array( self.get_fraction_exons_that_overlap_SNAP_from_QI() ) )
                              
    def get_AEDs_stats(self):
        return self.get_stats(np.array( self.get_AEDs() ) )
                              
    def get_eAEDs_stats(self):
        return self.get_stats(np.array( self.get_eAEDs() ) )


    def _parse_desc(self, desc):
        d = {}
        for e in desc:
            try:
                k,v = e.split('=')
                join = '='
            except: #Although not approp - I have seen some GFF desc terms sep by :
                k,v = e.split(':')
                join = ':'
            d[k] = v
        return d



class StringTieGFF(GFF):
    def __init__(self, allowable_sources=['StringTie'], allowable_types=['transcript', 'exon']):
        GFF.__init__(self, allowable_sources, allowable_types)
        self.exonidcnt = defaultdict(int)
    def add(self, line):
        ''' line = unsplit GFF line '''
        entry = line.strip().split('\t')
        if line:
            if line.startswith('#') and not self.headerpassed:
                self.header += line
            elif line.startswith('#') and self.headerpassed:
                pass #continue
            elif len(entry)>=9 and entry[1] in self.allowable_sources and entry[2] in self.allowable_types:
                # i.e. has descriptions column in col9
                desc = entry[8].rstrip(';').split(';')
                
                d = self._parse_desc(desc)
                if entry[2] == 'transcript': ## attributes are ID and geneID (analogous to Parent)
                    d['Parent'] = d['geneID']
                    if d['Parent'] not in self.genes.keys():
                        g = {}
                        g['ID'] = d['Parent']
                        g['Name'] = d['Parent']
                        gene = Gene(entry, attributes=g)
                        self.genes[g['ID']] = gene
                        self.num_genes += 1
                        #self.gene_lengths.append ( gene.length() ) ## For StringTie, this needs to done last
                    RNA = mRNA(entry, attributes=d)
                    self.num_mRNA += 1
                    self.mRNA_lengths.append( RNA.length() )
                    for gene in RNA.Parents():
                        self.genes[gene].add_mRNA(RNA)
                        self.RNA_parents[RNA.ID()].add( gene )
                elif entry[2] == 'exon': ## Attributes are Parent -- Parent corresponds to transcript
                    self.exonidcnt[d['Parent']] += 1
                    d['ID'] = d['Parent']+':exon_'+str(self.exonidcnt[d['Parent']])
                    exon = GenePart(entry, attributes=d)
                    self.num_exons += 1
                    self.exon_lengths.append( exon.length() )
                    genes = set([])
                    for RNA in exon.Parents():
                        for gene in self.RNA_parents[RNA]:
                            genes.add(gene)
                    for gene in list(genes): #should only be 1
                        self.genes[gene].add_exon( exon )



def write(f, l):
    f.write('\n'.join([str(e) for e in l]) + '\n')


def read_gff(args):
    if args.maker:
        gff = GFF()
        fh = args.maker
    elif args.stringtie:
        gff = StringTieGFF()
        fh = args.stringtie
    with open(fh) as gff_fh:
        for line in gff_fh:
            gff.add(line)
    return gff

def output(args, gff):    
    print "Num Genes", gff.get_num_genes()
    print "Num mRNA", gff.get_num_mRNA()
    print "Num Exons", gff.get_num_exons()
    print "Num Introns", gff.get_num_introns()
    print "ExonPerGeneStats", gff.get_exon_per_gene_stats()
    print "ExonLengthStats", gff.get_exon_length_stats()
    print "IntronLengthStats", gff.get_intron_length_stats()
    print "mRNALengthStats", gff.get_mRNA_length_stats()
    print "GeneLengthStats", gff.get_gene_length_stats()
    print "mRNA:Gene Ratio", gff.get_num_mRNA() / float( gff.get_num_genes() )
    print "Num mRNA perGeneStats (isoforms per gene)", gff.get_mRNA_per_gene_stats()
    if args.maker or args.full:
        print "Num CDS", gff.get_num_cds()
        print "Num 5pUTR", gff.get_num_5pUTR()
        print "Num 3pUTR", gff.get_num_3pUTR()
        print "Num genes with 5pUTR", gff.get_num_genes_with_5pUTR()
        print "Num genes with 3pUTR", gff.get_num_genes_with_3pUTR()
        print "CDSLengthStats", gff.get_cds_length_stats()
        print "5pUTRLengthStats", gff.get_5pUTR_length_stats()
        print "3pUTRLengthStats", gff.get_3pUTR_length_stats()
        print "Num5pUTRPerGeneStats", gff.get_num_5pUTR_per_gene_stats()
        print "Num3pUTRPerGeneStats", gff.get_num_3pUTR_per_gene_stats()
        print "Num Genes with Ontology Term", gff.get_num_genes_with_ontology()
        print "Num Genes with Blastp", gff.get_num_genes_with_blastp()
        print "Num Genes with Pfam", gff.get_num_genes_with_pfam()
        print "Num Genes with All 3 (intersect)", gff.get_num_genes_with_ontology_blastp_and_pfam()
        print "Num Genes with at least 1 of 3 (union)", gff.get_num_genes_with_ontology_blastp_or_pfam()
    if args.maker:
        print "Fraction Splice Sites Confirms By EST", gff.get_fraction_splice_sites_confirmed_by_EST_from_QI_stats()
        print "Fraction Splice Sites Confirms By Ab Initio", gff.get_fraction_splice_sites_confirmed_by_SNAP_from_QI_stats()
        print "Fraction Exons Overlap EST", gff.get_fraction_exons_that_overlap_EST_from_QI_stats()
        print "Fraction Exons Overlap EST or Protein", gff.get_fraction_exons_that_overlap_EST_or_protaln_from_QI_stats()
        print "Fraction Exons Overlap Ab Initio", gff.get_fraction_exons_that_overlap_SNAP_from_QI_stats()
        print "AED stats", gff.get_AEDs_stats()
        print "eAED stats", gff.get_eAEDs_stats() 


    if args.prefix:
        with open(args.prefix+'-Number-of-exons-per-gene.txt','w') as f:
            write(f, gff.get_exons_per_gene())
        with open(args.prefix+'-Number-of-mRNA-per-gene.isoformsPerGene.txt','w') as f:
            write(f, gff.get_mRNA_per_gene())
        with open(args.prefix+'-Gene-lengths.txt','w') as f:
            write(f, gff.get_gene_lengths())
        with open(args.prefix+'-pre-mRNA-lengths.txt','w') as f:
            write(f, gff.get_mRNA_lengths())
        with open(args.prefix+'-Exon-lengths.txt','w') as f:
            write(f, gff.get_exon_lengths())
        with open(args.prefix+'-Intron-lengths.txt','w') as f:
            write(f, gff.get_intron_lengths())

        if args.maker or args.full:
            with open(args.prefix+'-5pUTR-lengths.txt','w') as f:
                write(f, gff.get_5pUTR_lengths())
            with open(args.prefix+'-3pUTR-lengths.txt','w') as f:
                write(f, gff.get_3pUTR_lengths())
        if args.maker:
            with open(args.prefix+'-fraction_splice_sites_confirmed_by_EST.txt','w') as f:
                write(f, gff.get_fraction_splice_sites_confirmed_by_EST_from_QI())
            with open(args.prefix+'-fraction_splice_sites_confirmed_by_abinitio.txt','w') as f:
                write(f, gff.get_fraction_splice_sites_confirmed_by_SNAP_from_QI())
            with open(args.prefix+'-fraction_exons_that_overlap_EST.txt','w') as f:
                write(f, gff.get_fraction_exons_that_overlap_EST_from_QI())
            with open(args.prefix+'-fraction_exons_that_overlap_EST_or_protaln.txt','w') as f:
                write(f, gff.get_fraction_exons_that_overlap_EST_or_protaln_from_QI())
            with open(args.prefix+'-fraction_exons_that_overlap_abinitio.txt','w') as f:
                write(f, gff.get_fraction_exons_that_overlap_SNAP_from_QI())
            with open(args.prefix+'-AED-dist.txt','w') as f:
                write(f, gff.get_AEDs())
            with open(args.prefix+'-eAED-dist.txt','w') as f:
                write(f, gff.get_eAEDs()) 

    return

def run(args):
    output(args, read_gff(args))


## EXECUTE
run(args)

