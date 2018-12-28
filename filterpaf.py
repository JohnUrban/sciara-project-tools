#!/usr/bin/env python2.7
import sys
import argparse
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

 Return PAF lines given rules.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--paf', '-i',
                   type=str, required=True,
                   help='''Path to input PAF file.''')

parser.add_argument('--same', '-S',
                   action='store_true', default=False,
                   help='''Return lines where there was a "self-self" alignment as defined by some rules.
                        NOTE: this does not look at sequence names for the self-self definition.
                        It looks at sequence lengths and alignments.''')


parser.add_argument('--same_rules', '-SR',
                   type=str, default='1,1,1,1,1',
                   help='''Define "self-self" alignment with these rules.
                        Provide comma-sep list of: minLenRatio, maxLenRatio, minPctQueryMatch, minPctTargetMatch, minPctAlnMatch.
                        Default: 1,1,1,1,1.
                        Useful: 0.75,1.33,0.75,0.75,0.5.
                        It is treated as: (LenRatio within minLenRatio-maxLenRatio) OR (minPctQueryMatch AND minPctTargetMatch AND minPctAlnMatch)''')


parser.add_argument('--bubble', '-B',
                   action='store_true', default=False,
                   help='''Return lines where there the query was a "bubble" alignment as defined by some rules.
                        ''')


parser.add_argument('--bubble_rules', '-BR',
                   type=str, default='0.8,0.8,1',
                   help='''Define "bubble" alignment with these rules.
                        Provide comma-sep list of: minPctQueryInTarget,minPctAlnMatch,maxPctOfTarget.
                        Default: 0.8,0.8,1.
                        That is: at least 80% of query must be aligned in target with at least 80% matches in the alignment, and the query can be no longer than the target.
                        ''')

args = parser.parse_args()


class Paf(obj):
    def __init__(self, paffile):
        self.fmt = {0:str, 1:int, 2:int, 3:int, 4:str, 5:str, 6:int, 7:int, 8:int, 9:int, 10:int, 11:int}
        self.file = paffile
        self.key = {'query':0, 'qlen':1,'qstart':2,'qend':3,'strand':4,'target':5,'tlen':6,'tstart':7,'tend':8,'match':9,'alnlen':10,'mapq':11}
        with open(self.file) as fh:
            # only grabs columns 1-12
            self.paf = [[self.fmt[i](linelist[i]) for i in range(len(linelist))] for linelist in [line.strip().split()[:12] for line in fh.readlines()]]
        self.pafdf = pandas.DataFrame(paf, columns=['query', 'qlen','qstart','qend','strand','target','tlen','tstart','tend','match','alnlen','mapq']).sort_values(by='target','tstart','tend')
        self.iterpaf = iter(self.paf)

    def __iter__(self):
        return self

    def next(self):
        return self.iterpaf.next()

    def line2txt(self,line):
        return '\t'.join([str(e) for e in line])
    
paf = Paf(args.paf)

bubblerules = dict(zip(['minPctQueryInTarget','minPctAlnMatch','maxPctOfTarget'], [float(e) for e in args.bubble_rules.strip().split(',')]))

if args.same:
    rules = dict(zip(['minLenRatio', 'maxLenRatio', 'minPctQueryMatch', 'minPctTargetMatch', 'minPctAlnMatch'], [float(e) for e in args.same_rules.strip().split(',')]))
    for line in paf:
        qlen = float(line[paf.key['qlen']])
        tlen = float(line[paf.key['tlen']])
        match = float(line[paf.key['match']])
        alnlen = float(line[paf.key['alnlen']])
        pass_len_ratio = (tlen/qlen >= rules['minLenRatio'] or tlen/qlen <= rules['maxLenRatio'])
        pass_query_aln = (match/qlen >= rules['minPctQueryMatch'])
        pass_target_aln = (match/tlen >= rules['minPctTargetMatch'])
        pass_aln = (match/alnlen >= rules['minPctAlnMatch'])
        if pass_len_ratio or (pass_query_aln and pass_target_aln and pass_aln):
            print paf.line2txt(line)

