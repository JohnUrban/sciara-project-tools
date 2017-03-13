#!/usr/local/bin/python
import sys

if len(sys.argv) == 1:
    print "Usage: python scriptname ng/ul avgFragLen"
    quit()

concentration = sys.argv[1]
fragLength = sys.argv[2]

print float(concentration)*(1e+6)*(1/660.0)*(1/float(fragLength))
