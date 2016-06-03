from matplotlib import pyplot as plt

def gc(x):
    c=0
    for b in x:
        c += b.upper() in 'GC'
    return 100.0*c/len(x)

def complement(seq):
    ''' assumes uppercase '''
    c=""
    for b in seq:
            if b == "A": c += "T"
            elif b == "C": c += "G"
            elif b == "G": c += "C"
            elif b == "T": c += "A"
    return c    

def revcomp(seq):
    return complement(seq)[-1::-1]

def case_counter(seq):
    return sum(1 for b in seq if b.islower())



