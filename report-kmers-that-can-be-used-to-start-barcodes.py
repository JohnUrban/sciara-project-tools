from itertools import product

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


a = open("plasmid+A+B.seq",'r')
c = open("plasmid+C.seq", 'r')

A = ''
C = ''


for line in a:
    line = line.strip().split()
    for e in line:
            A+=e

for line in c:
    line = line.strip().split()
    for e in line:
            C+=e
a.close()
c.close()
A = A.upper()
C = C.upper()

six = get_all_kmers(6)

forA = [k for k in six if k not in A and rc(k) not in A and 'A' not in k and k[0] == "T"]
forC = [k for k in six if k not in C and rc(k) not in C and 'A' not in k and k[0] == "T"]


print "There are %d 6mers that do not appear (nor do their revcomps) in plasmid+A+B that do not have A and that start with T" % (len(forA))
print ("\n").join(forA)
print

print "There are %d 6mers that do not appear (nor do their revcomps) in plasmid+C that do not have A and that start with T" % (len(forC))
print ("\n").join(forC)
print


