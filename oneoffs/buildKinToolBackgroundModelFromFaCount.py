#!/usr/bin/env python2.7
import sys, argparse
import pandas as pd
import itertools as itr
from collections import defaultdict

parser = argparse.ArgumentParser(description="""
    
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('facount', metavar='facount', nargs='+',
                   type= str, 
                   help='''Path to facount.txt file (output of faCount or faCount.py).''')
parser.add_argument('--context', '-C', type=str, default="CG", help='''Middle Base to build a model for... Default CG (for now b/c that is why I developed it). Typically specified as A or C.''')
parser.add_argument('--klist', '-k', type=str, default='3,4,5,6,7', help='''Comma-separated list of kmer-sizes (integers) to build model for. Default: 3,4,5,6,7 ''')

args = parser.parse_args()





class FaCountMidBaseModel(object):
    def __init__(self, fh, kmers=[1,2,3,4,5,6,7], context='CG'):
        self.table = pd.read_csv(fh, sep="\t")
        self.kmers = kmers
        self.kdict = {k:[''.join(e) for e in itr.product('ACGT',repeat=k)] for k in self.kmers}
        self.model = {k:{i:defaultdict(int) for i in range(k)} for k in kmers}
        self.context = context
        self.clen = len(self.context)
        self.process_table()

    def process_table(self):
        for k in sorted(self.kdict.keys()):
            if k > self.clen: ## can't get context around anything shorter
                for kmer in self.kdict[k]:
                    for i in range(k-self.clen+1):
                        if kmer[i:i+self.clen] == self.context:
                            self.model[k][i][kmer] = self.table[kmer][0]
                
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
    


def run(args):
    print FaCountMidBaseModel(fh=args.facount[0], kmers=[int(e) for e in args.klist.split(',')], context=args.context)

try:
    run(args)
except IOError:
    pass
