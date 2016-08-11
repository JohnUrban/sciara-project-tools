#!/usr/bin/env python
import os, sys, argparse

parser = argparse.ArgumentParser(description="""

TODO

Deal with HPCdaligner with SLURM.

Assumes following in environment:
fasta2DB
DBsplit
HPCdaligner
daligner
LAsort
LAmerge

Assumes you are SLURM system had have following in environment:
sbatch

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--dbname', '-n',
                   type=str, required=True,
                   help='''Provide name for fasta2DB database that will be created and/or used.''')



args = parser.parse_args()

