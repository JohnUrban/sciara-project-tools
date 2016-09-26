#!/usr/local/bin/python
from collections import defaultdict, deque
import random, sys

## provide sampleSize and fragmentsFile
## sampleSize <= # fragments in fragments file
## Made 2013-2014


def sample(sampleSize, fragmentsFile, numFrags = 0):
    if numFrags == 0:
        f = open(fragmentsFile, 'r')
        for line in f:
            numFrags += 1
        f.close()
    sampledLineIndexes = random.sample(range(1,numFrags+1), sampleSize)
    sampledLineIndexes .sort()
    return deque(sampledLineIndexes )

def sampleBed(sampledLineIndexes, fragmentsFile):
    '''Takes line indexes sampled by sample() and the fragments file
    returns only the lines samples'''
    lineNum = 1
    f = open(fragmentsFile, 'r')
    while sampledLineIndexes:
        while lineNum != sampledLineIndexes[0]:
            lineNum += 1
            f.next()
        entry = f.next()
        lineNum += 1
        sys.stdout.write(entry)
        sampledLineIndexes.popleft()
    
    

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print "usage: *.py sampleSize fragmentsFile"
        quit()
    sampleSize = int(sys.argv[1])
    fragmestFile = sys.stdin
    sampledLineIndexes = sample(sampleSize, fragmentsFile)
    sampleBed(sampledLineIndexes, fragmentsFile)



    
