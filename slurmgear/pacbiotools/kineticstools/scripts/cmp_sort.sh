#!/bin/bash

export PATH=/users/jurban/software/localpy/pbh5tools/bin/:${PATH}

date

cmph5tools.py sort --deep ${IN_CMP} --tmpDir ${TMP_SORT}

date
