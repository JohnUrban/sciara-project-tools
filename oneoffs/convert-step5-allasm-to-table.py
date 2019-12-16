#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

DESCRIPTION -

    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--file', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file.''')

parser.add_argument('--rename_by_size', '-r', default=False, action='store_true', help=''' Return modified contig names after sorted by size longest to shortest.''')
parser.add_argument('--sort_by_size', '-s', default=False, action='store_true', help=''' Return contig info in order of size longest to shortest.''')


args = parser.parse_args()


d = dict()
size =defaultdict(list)
with open(args.file) as f:
    for line in f:
        #print line
        name, length = line.strip().split()
        #print name
        contents = name.split('_')
        #print contents
        contam = contents[0]
        num = contents[2]
        orig = "contig_"+num
        purgeclass = contents[-1]
        if 'contain' in name.lower():
            contain = 'yes'
            isrmtig='no'
        elif 'rmtig' in name.lower():
            contain='yes'
            isrmtig='yes'
        else:
            contain = 'no'
            isrmtig='no'
        redundant = 'yes'if 'redundant' in name.lower() else 'no'
            
        #print '\t'.join([orig, length, contam, isrmtig, contain, redundant, purgeclass])
        d[int(num)] = [orig, length, contam, isrmtig, contain, redundant, purgeclass, name]
        size[int(length)].append([orig, length, contam, isrmtig, contain, redundant, purgeclass, name])

if args.rename_by_size:
    print '#' + '\t'.join(['contig','length','tax_label', 'is_remove_tig', 'contains_rmtig', 'is_redundant', 'purge_haplotigs_class', 'orig_name', 'orig_info_name'])
    i=0
    for key in sorted(size.keys(), reverse=True):
        for contents in size[key]:
            i+=1
            #contents = size[key]
            print '\t'.join(['contig_'+str(i)]+contents[1:-1]+[contents[0]]+[contents[-1]])
elif args.sort_by_size:
    print '#' + '\t'.join(['contig','length','tax_label', 'is_remove_tig', 'contains_rmtig', 'is_redundant', 'purge_haplotigs_class', 'info_name'])
    for key in sorted(size.keys(), reverse=True):
        for contents in size[key]:
            print '\t'.join(contents)
else:
    print '#' + '\t'.join(['contig','length','tax_label', 'is_remove_tig', 'contains_rmtig', 'is_redundant', 'purge_haplotigs_class', 'info_name'])
    for key in sorted(d.keys()):
        print '\t'.join(d[key])
