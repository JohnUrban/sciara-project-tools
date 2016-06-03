#!/usr/bin/env python

import sys, random, argparse


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    N random barcodes of length L

    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument("-A", "--alphabet",
                   type=str, default='A,C,G,T',
                   help='''Comma-separated letters to make barcode with. Default: A,C,G,T''')
parser.add_argument("-L", "--length",
                   type=int, default=False,
                   help='''Barcodes should be length L''')
parser.add_argument("-N", "--number",
                   type=int, default=False,
                   help='''make N barcodes''')
parser.add_argument("-P", "--pairs",
                   action="store_true", default=False,
                   help='''make N barcode pairs''')
parser.add_argument("-w", "--nonself",
                   type=int, default=3,
                   help='''Nonself weight -- given last letter was b, higher weights give lower probablities that next letter will be b. Default: 2. w=1 gives uniform. w=0 gives homopolymers of length L.''')

parser.add_argument("-W", "--weights",
                    type=str, default=False,
                    help=''' Provide comma-separated list of integer weights for each letter in alphabet. If not used, they are weighted uniformly.''')

parser.add_argument("-C", "--cbind",
                   type=str, default=False,
                   help='''Combine seq output with columns from specified file.''')

parser.add_argument("-pre1", "--pre1",
                   type=str, default="",
                   help='''Add this prefix to all sequences generated. This is added to first in pair if --pairs specified.''')

parser.add_argument("-suf1", "--suf1",
                   type=str, default="",
                   help='''Add this suffix to all sequences generated. This is added to first in pair if --pairs specified.''')

parser.add_argument("-pre2", "--pre2",
                   type=str, default="",
                   help='''Add this prefix to all sequences generated. This is added to second in pair if --pairs specified.''')

parser.add_argument("-suf2", "--suf2",
                   type=str, default="",
                   help='''Add this suffix to all sequences generated. This is added to second in pair if --pairs specified.''')

parser.add_argument("-gc", "--gc",
                   type=str, default="0,100",
                   help='''Only return sequences within this range of GC content. Default: 0,100. Provide 2 comma-sep integers. e.g. 60,67.
Note: as it is now, this means you may not get back N barcodes. You will get back the subset of N generated that fits into the GC constraints.''')


args = parser.parse_args()




alphabet = args.alphabet.strip().split(",")

weights = {k:1 for k in alphabet}

if args.weights:
    args.weights = [int(e) for e in args.weights.strip().split(",")]
    for i in range(len(alphabet)):
        weights[alphabet[i]] = args.weights[i]


def getmodel(b,nonself_weight,alphabet, A_weight=1,C_weight=1,G_weight=1,T_weight=1):
    if b in "aA":
        return A_weight*"A" + C_weight*"C"*nonself_weight + G_weight*"G"*nonself_weight + T_weight*"T"*nonself_weight
    elif b in "cC":
        return A_weight*"A"*nonself_weight + C_weight*"C" + G_weight*"G"*nonself_weight + T_weight*"T"*nonself_weight
    elif b in "gG":
        return A_weight*"A"*nonself_weight + C_weight*"C"*nonself_weight + G_weight*"G" + T_weight*"T"*nonself_weight
    elif b in "tT":
        return A_weight*"A"*nonself_weight + C_weight*"C"*nonself_weight + G_weight*"G"*nonself_weight + T_weight*"T"

def getmodel(b,nonself_weight,alphabet):
    return b + ('').join([e*nonself_weight for e in alphabet if e != b])

def getmodel(b,nonself_weight,alphabet, weights):
    return b*weights[b] + ('').join([e*nonself_weight*weights[e] for e in alphabet if e != b])

def gc(x):
    c=0
    for b in x:
        if b in 'GC':
            c+=1
    return 100.0*c/len(x)

if args.cbind:
    lines = []
    for line in open(args.cbind):
        lines.append(line.strip())
    if lines[0][0] == "#":
        print lines[0]
        lines.pop(0)
    args.number = len(lines)

def comp(s):
    n = ''
    for b in s:
        if b in 'A':
            n+='T'
        elif b in 'C':
            n+='G'
        elif b in 'G':
            n+='C'
        elif b in 'T':
            n+='A'
    return n

def rc(s):
    return comp(s)[-1::-1]

def get_all_kmers(k):
    assert k > 0 and k < 7
    if k == 1:
        return [''.join(e) for e in product("ACGT")]
    elif k == 2:
        return [''.join(e) for e in product("ACGT","ACGT")]
    elif k == 3:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT")]
    elif k == 4:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT","ACGT")]
    elif k == 5:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT","ACGT","ACGT")]
    elif k == 6:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT","ACGT","ACGT","ACGT")]
    
init = alphabet
gclimit = [int(e) for e in args.gc.strip().split(",")]

for i in range(args.number):
    seq = init[random.randint(0,len(init)-1)]
    model = getmodel(seq[-1],args.nonself,alphabet,weights)
    for j in range(1,args.length):
        seq += model[random.randint(0,len(model)-1)]
        model = getmodel(seq[-1],args.nonself,alphabet,weights)
    seqgc = gc(args.pre1+seq+args.suf1)
    if args.pairs:
        seq2 = init[random.randint(0,len(init)-1)]
        model = getmodel(seq[-1],args.nonself,alphabet,weights)
        for j in range(1,args.length):
            seq2 += model[random.randint(0,len(model)-1)]
            model = getmodel(seq2[-1],args.nonself,alphabet,weights)
        seqgc2 = gc(args.pre2+seq2+args.suf2)
        if seqgc >= gclimit[0] and seqgc <= gclimit[1] and seqgc2 >= gclimit[0] and seqgc2 <= gclimit[1]:
            if args.cbind:
                print lines[i] + "\t" + args.pre1+seq+args.suf1 + "\t" + args.pre2+seq2+args.suf2
            else:
                print args.pre1+seq+args.suf1 + "\t" + args.pre2+seq2+args.suf2
    else:
        if seqgc >= gclimit[0] and seqgc <= gclimit[1]:
            if args.cbind:
                print lines[i] + "\t" + args.pre1+seq+args.suf1
            else:
                print args.pre1+seq+args.suf1
