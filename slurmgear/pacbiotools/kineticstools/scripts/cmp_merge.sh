#!/bin/bash

export PATH=/users/jurban/software/localpy/pbh5tools/bin/:${PATH}

date

cmph5tools.py merge --outFile ${MERGEDCMP} ${CMPDIR}/*.cmp.h5

date
