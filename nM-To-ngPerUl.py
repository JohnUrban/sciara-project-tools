#!/usr/local/bin/python
import sys

if len(sys.argv) == 1:
    print "Usage: python scriptname nM avgFragLen"
    quit()

nM = sys.argv[1]
fragLength = sys.argv[2]

print float(nM)*(1e-6)*(660.0)*(float(fragLength))
