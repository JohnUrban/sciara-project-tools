#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict, Counter
from Bio import SeqIO
from numpy.random import choice


parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Take in a FASTA, learn a background model for kineticsTools.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('fasta', metavar='fasta', nargs='*',
                   type= str, 
                   help='''Path to as many FASTAs as need be.''')


parser.add_argument('--midbase', '-b', type=str, default="ACGT", help='''Middle Base to build a model for... Default (all): ACGT. Typically specified as A or C.''')
parser.add_argument('--windowsize', '-w', type=int, default=41, help='''Scan sequences in FASTA with this window size. Default: 41 (the seqlen reported in kintools outputs).''')
parser.add_argument('--mididx', '-i', type=int, default=20, help='''Middle Base Index. Default = 20. Usually don't need to change this unless window size is changed.''')
parser.add_argument('--klist', '-k', type=str, default='1,2,3,4,5,6,7', help='''Comma-separated list of kmer-sizes (integers) to build model for. Default: 1,2,3,4,5,6,7 ''')
parser.add_argument('--scores', '-s', type=str, default=False, help='''Provide single-column text file of scores. It will be read in. A score will be randomly sampled from the list for each window. The score will be used to create a randomly weighted background model.''')
parser.add_argument('--norevcomp', '-N', action='store_true', default=False, help='''Do not look at reverse complements of given sequences.''')

args = parser.parse_args()






class MidBaseModel(object):
    def __init__(self, kmers=[1,2,3,4,5,6,7], mid=20):
        self.kmers = kmers
        self.model = {k:{i:defaultdict(int) for i in range(k)} for k in kmers}
        self.mid = mid

    def add(self, seq, score=1):
        '''seq is string of upper-case'''
        seq = seq.upper()
        for k in self.kmers:
            for i in sorted(self.model[k].keys()):
                left = self.mid - i
                right = self.mid + (k-i)
                self.model[k][i][seq[left:right]] += score

    def __str__(self):
        #return str(self.model)
        out = ''
        gate = set('ACGT') ## does not allow kmers with N or any other non-canonical
        for k in self.kmers:
            for i in sorted(self.model[k].keys()):
                for seq in sorted(self.model[k][i].keys()):
                    chars = set(seq)
                    if not chars.difference(gate):
                        ## Add 1 to i so it is 1-based
                        out += '\t'.join([str(e) for e in [k, i+1, seq, self.model[k][i][seq]]]) + '\n'
        return out.strip()
    
#### EXECUTE







try:
    ## Initialize model
    midbasemodel = MidBaseModel(kmers=[int(e) for e in args.klist.split(',')], mid=args.mididx)
    
    ## Initialize scoring 
    score = 1 ## doesn't change unless args.scores used
    if args.scores:
        scores, counts = [np.array(tup) for tup in zip(*Counter([int(e.strip()) for e in args.scores.readlines()]))]
        nscores = float(counts.sum())
        probs = counts/nscores

    ## Loop over fastas
    for fh in args.fasta:

        ## Loop over sequences in fasta
        for fa in SeqIO.parse(fh, 'fasta'):

            ## Loop over windows in sequence
            for pos in range(len(fa)-args.windowsize+1):
                
                ## Get subset window
                subfa = fa[pos:pos+args.windowsize]

                ## Add fwd seq if correct midbase
                seq = str(subfa.seq)
                if seq[args.mididx] in args.midbase:
                    if args.scores:
                        score = choice(a = scores, size = 1, replace = True, p = probs)
                    midbasemodel.add(seq, score=score)

                ## Add RevComp seq
                if not args.norevcomp:
                    seq = str(subfa.seq.reverse_complement())
                    if seq[args.mididx] in args.midbase:
                        if args.scores:
                            score = choice(a = scores, size = 1, replace = True, p = probs)
                        midbasemodel.add(seq, score=score)

except IOError:
    pass



## RETURN MODEL
print midbasemodel


quit()

