from collections import defaultdict

enzymes = {'ApaI':'GGGCCC','BamHI':'GGATCC','EcoRI':'GAATTC','EcoRV':'GATATC','HindIII':'AAGCTT', 'KpnI':'GGTACC', 'PstI':'CTGCAG', 'XhoI':'CTCGAG', 'SalI':'GTCGAC', 'XbaI':'TCTAGA', 'BglI':'GCCNNNNNGGC', 'MstII':'CCTNAGG', 'HaeII':'RGCGCY'}
enzymes['Sau3AI'] = 'GATC'
enzymes['HpaII'] = 'CCGG'
enzymes['MspI'] = 'CCGG'

## I need this here until I implement regular expressions
enzymes['BglI_1'] = 'GCCAATTCGGC'
enzymes['BglI_2'] = 'GCCCAGTGGGC' #GCCC AGC TGGC
enzymes['BglI_3'] = 'GCCCAGCTGGC'

enzymes['BssSI'] = 'CACGAG'
enzymes['BssSIrc'] = 'CTCGTG'

enzymes['BspQI'] = 'GCTCTTC'
enzymes['BspQIrc'] = 'GAAGAGC'


## These below are needed for output stage (BED)
enzymes['5p'] = '*'
enzymes['3p'] = '*'

enzymelookup = defaultdict(list) ## THIS IS FILLED OUT ON THE FLY - I've added this in anticipation of the Regex approach of searching a sequence in one pass


## Dinucleotides
enzymes['CG'] = 'CG'
enzymes['AG'] = 'AG'
enzymes['AGrc'] = 'CT'

## Trinucleotides
enzymes['GCG'] = 'GCG'
enzymes['GCGrc'] = 'CGC'
enzymes['GAG'] = 'GAG'
enzymes['GAGrc'] = 'CTC'

