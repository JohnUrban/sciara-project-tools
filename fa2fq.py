#!/usr/bin/env python
## June 9, 2014
import sys

def fastaToFastq(faFileConnection):
    ''' fastq structure is 4 lines per unit
        lines 1-4 == @name, sequence, +name, quality scores
        This extracts just the sequences (line 2s)
        it prints each so as to not store in memory
        Thus, this needs to be used from command line and redirected
        --
        EXPECTS FA file to be in 2 line per entry format (name:seq pairs)
        Use fasta_formatter from fastX toolkit to re-format fa file if nec.'''
    ctr = 0
    pairCount = 0
    for line in faFileConnection:
        ctr += 1
        pairCount += 1
        if ctr%2 == 1:
            print "@"+ line[1:-1]
        if ctr%2 == 0:
            seqLen = len(line[:-1])
            print line[:-1]
        if pairCount == 2:
            print "+"
            print 'J'*seqLen
            pairCount = 0


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print
        print "Use 1: fa2fq.py file.fa"
        print "Use 2: other commands | fa2fq.py"
        print "Use 3: other commands | fa2fq.py -"
        print
    elif sys.argv[1] == "stdin" or sys.argv[1] == "-":
        fastaToFastq(sys.stdin)
    elif len(sys.argv) == 2:
        connection = open(sys.argv[1], 'r')
        fastaToFastq(connection)
