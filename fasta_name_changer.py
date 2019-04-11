#!/usr/bin/env python

import sys, argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description="""

DESCRIPTION - version 24Mar2015
    
    Change/modify name of each fasta entry in a file with MANY fasta entries.

    Fasta names usually follow ">name" format.

    Some may follow ">name otherinfo".
    This can screw up some programs.
    A way to keep only the part up to the first white space is:
       $ fasta_name_changer.py -f file.fa -k 1
    Or keeping only the second block of info as the name:
       $ fasta_name_changer.py -f file.fa -k 2
    etc.

    Add nameinfo in front of name or selected part of name:
       $ fasta_name_changer.py -f file.fa -F prefix
       $ fasta_name_changer.py -f file.fa -k 1 -F prefix

    Add nameinfo in  back of name or selected part of name:
       $ fasta_name_changer.py -f file.fa -B prefix
       $ fasta_name_changer.py -f file.fa -k 1 -B prefix

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str, default=False, required=True,
                   help='''Path to input file in FASTA format.''')
parser.add_argument('--keep', '-k',
                   type=int, default=False,
                   help='''Break name up by whitespace. Keep only kth chunk (e.g. 1).''')
parser.add_argument('--front', '-F',
                   type=str, default=False,
                   help='''Add this info to front of name (operation occurs after -k)''')
parser.add_argument('--back', '-B',
                   type=str, default=False,
                   help='''Add this info to back of name (operation occurs after -k)''')
parser.add_argument('--replace', '-r',
                   type=str, default=False,
                   help='''replace name with given string -- often used with --number to append entry number to replacement word.''')
parser.add_argument('--number', '-n',
                   action='store_true', default=False,
                   help='''Add entry number (e.g. _1) to end of name. Often sed with --replace.''')
parser.add_argument('--join', '-j',
                   action='store_true', default=False,
                   help='''Break name/description at white space and join with underscore.''')

parser.add_argument('--joinchar', '-J',
                   type=str, default='_',
                   help='''Used only with --join. Change default join character of underscore to something else. e.g -J , ''')

parser.add_argument('--key', default=False, type=str, help=''' Provide filename to store new_name-to-old_name map file (3 columns: new name, old name, old description).''')
parser.add_argument('--sort', default=False, action='store_true', help=''' Before renaming, sort sequences from shortest to longest. Also specify --revsort to do longest to shortest.''')
parser.add_argument('--namesort', default=False, action='store_true', help=''' Ue with --sort to do sort names instead.''')
parser.add_argument('--revsort', default=False, action='store_true', help=''' Ue with --sort to do longest to shortest or with --sort and --namesort to do reverse order of names.''')


args = parser.parse_args()
if args.keep:
    chunk = args.keep-1
if args.fasta == "-" or args.fasta == "stdin":
    args.fasta = sys.stdin
if args.key:
    key = open(args.key,'w')
## SeqIO automatically takes ">name" from ">name other info" with entry.name
## entry.description gives whole name "name other info" -- removes ">" in front



## Only load entire file into memory if you need to sort
if args.sort:
    seqlist = []
    for entry in SeqIO.parse(args.fasta, "fasta"):
        desc = str(entry.description)
        name = str(entry.name)
        seq = str(entry.seq)
        seqlen = len(seq)
        seqlist.append([desc, seq, seqlen, name])
    

    if args.sort:
        if args.namesort:
            sort_i = 0
        else:
            sort_i = 2
        seqlist.sort(key = lambda x: x[sort_i], reverse = args.revsort) ## sort on length

if args.sort:
    SS=seqlist
else:
    SS=SeqIO.parse(args.fasta, "fasta")


seq_n=0
for entry in SS:   
#for entry in seqlist:
#for entry in SeqIO.parse(args.fasta, "fasta"):
    seq_n += 1
    #name = entry.description
    #name = entry[0]
    name = entry[0] if args.sort else entry.description
    name = (args.joinchar).join(name.split()) if args.join else name
    if args.keep:
##        i = args.keep-1
        name = name.split()[chunk]
    if args.front:
        name = args.front+name
    if args.back:
        name = name+args.back
    if args.replace:
        name = args.replace
    if args.number:
        name += "_"+str(seq_n)
    #name = ">"+name
    sys.stdout.write(">" + name + "\n")
    outseq = entry[1] if args.sort else str(entry.seq)
    sys.stdout.write(outseq + "\n")
    if args.key:
        outkey = name+'\t'+entry[3]+'\t'+entry[0]+'\n' if args.sort else name+'\t'+entry.name+'\t'+entry.description+'\n'
        key.write(outkey)

if args.key:
    key.close()
