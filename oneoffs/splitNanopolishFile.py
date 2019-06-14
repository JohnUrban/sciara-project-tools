#!/usr/bin/env python2.7

import sys, argparse


parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    Split nanopolosh file up.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i', type=str, required=True, help='''Path to nanopolish file.''')

args = parser.parse_args()



    
#### EXECUTE

model = 'nevergoingtoseethis'
kmer = 'nevergoingtoseethis'
outfile = 'nergoingtoseethis'

def close(g):
    try:
        g.close()
    except:
        pass

try:
    ##
    with open(args.input) as f:
        for line in f:
            if line:
                if line.startswith('#') or line.startswith('model'):
                    continue
                splitline = line.strip().split()
                linemodel = splitline[0]
                linekmer = splitline[1]
                if linemodel != model or linekmer != kmer:
                    close(outfile)
                    model = linemodel
                    kmer = linekmer
                    newfile = model + '.' + kmer + '.txt'
                    outfile = open(newfile, 'w')
                    outfile.write(line.strip())
                else:
                    outfile.write(line.strip())
    close(outfile)
    

except IOError:
    pass ## Fail silently here (e.g. when piped into less or head)




quit()

