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
    return 100.0*sum( [1 for b in x if b in 'GCgc'] )/len(x)



def comp(seq):
    intab=  'actguACTGUNRY'
    outtab= 'tgacaTGACANYR'
    transtab = string.maketrans(intab, outtab)
    return seq.translate(transtab)

def revcomp(seq):
    return comp(seq)[-1::-1]



def get_all_kmers(k):
    return [''.join(e) for e in [e for e in product('ACGT', repeat=k)]]
    
