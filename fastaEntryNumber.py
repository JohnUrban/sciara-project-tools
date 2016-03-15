#!/usr/local/bin/python
import sys


if len(sys.argv) == 1:
    print
    print "Usage: fastaEntryNumber.py file.fa entryName"
    print "Assumes you gave exact/precise entry name to find"
    print "Returns the entrty number of only the first time it appears"
    print
    quit()


f = open(sys.argv[1], 'r')
entry = sys.argv[2]

count = 0
for line in f:
    if line[0] == ">":
        count += 1
        line = line[1:-1]
        if entry == line:
            print count
            quit()
