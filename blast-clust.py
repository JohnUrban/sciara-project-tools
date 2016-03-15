#!/usr/bin/env python
import sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="""
    By John Urban. Trying to recapitulate/inspired by: http://www.ncbi.nlm.nih.gov/Web/Newsltr/Spring04/blastlab.html
    (Could not find the code to blastclust though)
    Schatz uses on nanopore assembly paper -- -b F -p F -e F -L 0.80 -S 60 -W 14
    Using above link and http://etutorials.org/Misc/blast/Part+V+BLAST+Reference/Chapter+13.+NCBI-BLAST+Reference/13.9+blastclust+Parameters/
    Their parameters seem to mean:
    -L 0.8 -b F (80% of length of one contig or the other has to be covered --- -b F seems to mean "both=False")
    -p F -> protein=False --> nucleotide
    -S 60 -- requires 60% identity
    -W word_size=14 -- I supposed for initial blast seed
    -e F (default is False anyway) -- this is just something that enables ID parsing --doesnt matter
    So try blasting with -word_size=14, -qcov 80, -perc_identity 60
    
    Instructions:
    1. Make fasta file of sequences to be clustered.
    2. use makeblastdb to make blastdb of them: makeblastdb -in targets.fa -dbtype nucl -out targets
    3. use blastn with fasta file against its own database to obtain all blast hit information.
        use promiscuous/non-stringent parameters for:
            percent identity (e.g. 70%)
            query coverage (qcov_hsp_perc) (e.g. 70%)
        other things can be used to speed up search for large database:
            evalue can be 1e-50 or 1e-100 (1e-5 if promiscuous)
            word_size can be >= 15 to speed up search, especially when plan to require high percent identity for clustering
        use output format 6 and use following ordering and terms:
            '6 qseqid sseqid pident length qlen slen evalue sstrand qcovs qcovhsp qstart qend sstart send'
        example blastn:
        blastn -qcov_hsp_perc $QCOV -perc_identity $PID -task blastn -db $DB -query $QUERY -evalue $E -word_size $WORD -gapopen 2 -gapextend 2 -penalty -3 -reward 2 -dust no -outfmt '6 qseqid sseqid pident length qlen slen evalue sstrand qcovs qcovhsp qstart qend sstart send' > allhit.txt
    4. use blast-clust.py on 
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument("-b", "--blastfile",
                   type= str, default=False, required=True,
                   help='''Tab-delim blast output file.''')
parser.add_argument("-q", "--qcov",
                   type=float, default=100, 
                   help='''What percentage of the query should be covered by subject. Default: 100
 Note: often ends of the query will not be included in the alignment if they do not match the subject well and therefore qcov will not be 100%. Try 99% or 95% to allow this.''')
parser.add_argument("-p", "--pctid",
                   type=float, default=100, 
                   help='''What should the minimum percent identity be to be included in cluster. Default: 100''')

args = parser.parse_args()



class BlastEntry(object):
    def __init__(self, entry):
        #entry is simply a line from the blast file, unadulterated -- keep as string
        # 'qseqid sseqid pident length qlen slen evalue sstrand qcovs qcovhsp qstart qend sstart send'
        entry = entry.strip().split()
        self.entry = entry
        self.query = entry[0]
        self.subject = entry[1]
        self.pctid = float(entry[2])
        self.length = int(entry[3])
        self.qlen = int(entry[4])
        self.slen = int(entry[5])
        self.evalue = float(entry[6])
        self.strand = entry[7]
        self.qcov = float(entry[8])
        self.qcovhsp = float(entry[9])
        self.qstart = int(entry[10])
        self.qend = int(entry[11])
        self.sstart = int(entry[12])
        self.send = int(entry[13])

    def can_be_clustered(self, qcovcutoff=100, pctidcutoff=100):
        if self.subject != self.query:
            if self.slen >= self.qlen: ##only allow subjects to contain same_size or smaller queries in their clusters
                if self.qcov >= qcovcutoff and self.pctid >= pctidcutoff: ##not entirely sure if I should use qcov or qcovhsp for this cutoff
                    return True
        return False
        
blastout = open(args.blastfile, 'r')
clusters = defaultdict(list)
for entry in blastout:
    entry = BlastEntry(entry)
    if entry.can_be_clustered(qcovcutoff=args.qcov, pctidcutoff=args.pctid):
        clusters[entry.subject].append(entry.query)

for subject in clusters.keys():
    sys.stdout.write(subject + ":\t" + (", ").join(clusters[subject])+"\n")
    
    
