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
