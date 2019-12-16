###!/usr/bin/env python2.7
##
##import sys, argparse
##from collections import defaultdict, Counter
###from Bio import SeqIO
###from numpy.random import choice
##import scipy
##from scipy.stats import chisquare
##
##parser = argparse.ArgumentParser(description="""
##
##DESCRIPTION -
##
##    Read in a test result,
##    Output filtered lines.
##
##    """, formatter_class= argparse.RawTextHelpFormatter)
##
##parser.add_argument('--minK', '-K', type=int, default=2, help='''Minimum kmer length to output. Default = 2.''')
##signif  = parser.add_mutually_exclusive_group()
##signifparser.add_argument('--pvalcutoff', '-p', type=float, help='''P-value cutoff. Use 1 to return all.''')
##signifparser.add_argument('--qvalcutoff', '-q', type=float, help='''FDR cutoff. All p-values are transformed into q-values via BH procedure.''')
##residual = parser.add_mutually_exclusive_group()
##residual.add_argument('--minResidual', '-r', type=float, help='''Minimum residual cutoff. When this is specified, all residuals are treated the same regardless of test group. ''')
##residual.add_argument('--minMAD', '-R', type=float, help='''Minimum numer of MADs above Median. When this is specified, the distribution of residuals within a group is used to determine the enriched kmers within that group. ''')
##parser.add_argument('--klist', '-k', type=str, default='1,2,3,4,5,6,7', help='''Comma-separated list of kmer-sizes (integers) to build model for. Default: 1,2,3,4,5,6,7 ''')
##
##args = parser.parse_args()
##
##
##
##
##
##
##class KinTest(obeject):
##    def __init__(self, fh, kmers=[1,2,3,4,5,6,7]):
##        self.fh
##        self.model = {k:{i:defaultdict(int) for i in range(k)} for k in kmers}
##        self.__process()
##    def __process(self):
##        self.lines
##        with open(self.fh) as f:
##            
##
##    
###### EXECUTE
##
##try:
##    ## Initialize models
##    testing = CompareModels(test_fh = args.test, bg_fh = args.background, kmers=[int(e) for e in args.klist.split(',')], pseudo=args.pseudo)
##
##    ## RETURN MODEL
##    print testing
##
##except IOError:
##    pass ## Fail silently here (e.g. when piped into less or head)
##
##
##
##
##quit()
##
