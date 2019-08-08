#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="""

DESCRIPTION - clean up GFF "Parent" slot.


Use case:
When using Maker2 for annotation one can have it output all models: keep_pred=1.
Then one can use interproscan to find functional domains.
Then one can take all with AED < 1 and/or with a functional domain (standard set).

Problem:
I found that it doesn't fully clean up Parent issues.
There could be a gene with >1 isoform - say 2 for now.
Only 1 isoform is kept due to having a functional domain.
That isoform may have had exons for example that had more than 1 isoform as parents.
That parent slot is not cleaned -- and NCBI validation throws errors.

Example:
- Parent=maker-contig_12-snap-gene-3.702-mRNA-1,maker-contig_12-snap-gene-3.702-mRNA-2
- When only maker-contig_12-snap-gene-3.702-mRNA-1 was retained.


This script will read in all of the IDs found when column 2 = "maker" or given word.
It will then go through and print out the GFF again.
When it encounters column2=maker lines, it will only return parents that were in IDs.

(Exact matching)

Note for future:
- Since the ID should show up before it appears as a Parent, this can probably be done in a single pass.

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--gff', '-g', '-i', '-f',
                   type= str, default=False, required=True,
                   help='''Path to GFF.''')

parser.add_argument('--word', '-w', 
                   type=str, default='maker', 
                   help='''What word to filter on in column2. Set to ALLWORDS if you want to ignore columns to and include everything.''')

args = parser.parse_args()




def not_a_target_line(line):
    ''' line = a gff line stripped and split already.'''
    if args.word == 'ALLWORDS':
        return False #all are target lines
    return line[1] != args.word


def is_a_target_line(line):
    ''' line = a gff line stripped and split already.'''
    return (line[1] == args.word or args.word == 'ALLWORDS')

def assert_is_a_target_line(line):
    try:
        assert is_a_target_line(line)
    except AssertionError:
        print line
        print args.word, line[1], 'ALLWORDS'
        print "AssertionError: Column2 did not match -w word argument. See above (args.word, column2, ALLWORDS)."
        quit()
    
def return_clean_parent(k, join, v):
    alt_v = []
    for sub_v in v.split(','):
        if sub_v in IDs:
            alt_v.append(sub_v)

    new_v = ','.join(alt_v)
    return k + join + new_v

def process_target_line(line):
    # Split up description column.
    desc = line[8].rstrip(';').split(';')
    altdesc = []
    # Make kv dictionary
    d = {}
    for e in desc:
        try:
            k,v = e.split('=')
            join = '='
        except: #GO terms, etc sep by :
            k,v = e.split(':')
            join = ':'
        d[k.strip('_')] = v

        if k == 'Parent':
            altdesc.append( return_clean_parent(k, join, v) )
        else:
            altdesc.append( e )

    # Return to string formatting
    altdesc = ';'.join(altdesc)

    # Update line
    line[8] = altdesc
    return line
        

## COLLECT IDs
IDs = set([])
with open(args.gff) as f:
    for line in f:
        if line:
            line = line.strip().split('\t')
            if not_a_target_line(line):
                continue
            if line[0][0] != '#' and len(line)>3:
                desc = line[8].rstrip(';').split(';')
                for e in desc:
                    try:
                        k,v = e.split('=')
                    except: #GO terms, etc sep by :
                        k,v = e.split(':')
                    if k == 'ID':
                        IDs.add( v )
                        break
                    


                    
## RETURN WITH CLEANED PARENT SLOTS
with open(args.gff) as f:
    for line in f:
        if line:
            # Return commented lines (usually header) as is
            if line.startswith('#'):
                print line.strip()
                continue

            # Strip/Split line
            line = line.strip().split('\t')

            # Return non-target lines as is
            if not_a_target_line(line):
                print '\t'.join(line)
                continue

            # All else should be target lines here, but check anyway
            assert_is_a_target_line(line)

            # Process - the len(line)>3 check was propgated from previous scripts where I guess I found issues...
                # For now I will comment it out and assume GFFs should be correct
            #   if len(line)>3:

            # Process target line
            line = process_target_line(line)
            print '\t'.join(line)

