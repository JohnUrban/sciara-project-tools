#!/usr/bin/env python
import sys, random, argparse
import numpy as np
from makeRandomSeqTools import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION

    Don't use this.

    So many better tools!

    I just needed this for a moment.

        

    """,
    formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-O', '--refpre',
                    type=str, default='reference',
                    help='''Default reference.''')

parser.add_argument('-o', '--readpre',
                    type=str, default='reads',
                    help='''Default reads.''')

parser.add_argument('-L', '--reflen',
                    type=int, default=100000,
                    help='''Default 100000.''')
parser.add_argument('-N', '--nrefs',
                    type=int, default=1,
                    help='''Default 1.''')

parser.add_argument('-n', '--nreads',
                    type=int, default=10,
                    help='''Default 10.''')

parser.add_argument('-m', '--meanreadlen',
                    type=int, default=5000,
                    help='''Default 5000.''')

parser.add_argument('-s', '--sdreadlen',
                    type=int, default=500,
                    help='''Default 500.''')

parser.add_argument('-l', '--minreadlen',
                    type=int, default=500,
                    help='''Default 500.''')



parser.add_argument("-A", "--alphabet",
                   type=str, default='A,C,G,T',
                   help='''Comma-separated letters to make barcode with. Default: A,C,G,T''')
parser.add_argument("-w", "--nonself",
                   type=int, default=3,
                   help='''Nonself weight -- given last letter was b, higher weights give lower probablities that next letter will be b. Default: 2. w=1 gives uniform. w=0 gives homopolymers of length L.''')

parser.add_argument("-W", "--weights",
                    type=str, default=False,
                    help=''' Provide comma-separated list of integer weights for each letter in alphabet. If not used, they are weighted uniformly.''')

args = parser.parse_args()




alphabet = args.alphabet.strip().split(",")

weights = {k:1 for k in alphabet}

if args.weights:
    args.weights = [int(e) for e in args.weights.strip().split(",")]
    for i in range(len(alphabet)):
        weights[alphabet[i]] = args.weights[i]



init = alphabet

with open(args.refpre+'.fasta','w') as f:
    d={}
    for i in range(args.nrefs):
        f.write( '>'+args.refpre+'_'+str(i)+'\n')
        seq = init[random.randint(0,len(init)-1)]
        model = getmodel(seq[-1],args.nonself,alphabet,weights)
        for j in range(1,args.reflen):
            seq += model[random.randint(0,len(model)-1)]
            model = getmodel(seq[-1],args.nonself,alphabet,weights)
        f.write(seq+'\n')
        d[i] = seq


with open(args.readpre+'.fasta','w') as f:
    for i in range(args.nreads):
        rlen = int( np.random.normal(args.meanreadlen,args.sdreadlen) )
        while rlen < args.minreadlen:
            rlen = int( np.random.normal(args.meanreadlen,args.sdreadlen) )
        refsel = np.random.choice(d.keys())
        start = int( np.random.uniform(0, len(d[refsel])-rlen) )
        f.write('>'+args.readpre+'_'+str(i)+'\n'+d[refsel][start:start+rlen]+'\n')
        
        
