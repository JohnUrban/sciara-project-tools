#!/usr/local/bin/python
from collections import defaultdict
import sys

## Made 2013-2014
## provide genome file and motif centers file

def siteBoundaries(genomeFile, motifCentersFile):
    ''' provide genome file (chr names and lengths)
    and bed file of centers (chr start end) -- can be just chr start'''
    g = open(genomeFile, 'r')
    centers = open(motifCentersFile, 'r')
    sitebounds = dict()
    ## enter 0 and length for each chr into dict
    for entry in g:
        entry = entry[:-1].split()
        chrom = entry[0]
        length = int(entry[1])
        sitebounds[chrom] = [0, length]
    g.close()
    for entry in centers:
        entry = entry[:-1].split()
        chrom = entry[0]
        start = int(entry[1])
        sitebounds[chrom].append(start)
    centers.close()
    #sort
    for chrom in sitebounds.keys():
        sitebounds[chrom].sort()
    return sitebounds

def fragmentsBed(sitebounds):
    '''sitebounds is utput dictionary from siteBoundaries()
    this function write out bed file of all fragments defined therein
    A fragment is the space between site boundaries -- e.g. 0 to first non-zero number'''
    try:
        for chrom in sitebounds.keys():
            for startIndex in range(len(sitebounds[chrom])-1):
                bedEntry = chrom + "\t" + str(sitebounds[chrom][startIndex]) + "\t" + str(sitebounds[chrom][startIndex+1]) + "\n"
                sys.stdout.write(bedEntry)
    except IOError:
        pass
            
    
    

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print "usage: *.py genomeFile motifCentersFile"
        quit()
    genomeFile = sys.argv[1]
    motifCentersFile = sys.argv[2]
    sitebounds = siteBoundaries(genomeFile, motifCentersFile)
    fragmentsBed(sitebounds)
