#!/usr/bin/env python2.7
import sys, os, datetime
import argparse
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from random import choice
from time import sleep

parser = argparse.ArgumentParser(description="""

Prior to this:
NT=~/software/blast/db/nt
NR=~/software/blast/db/nr
REP=../../../repeatmodeler/sciara_canu-families.fa
GFF=~/scratch/repeatmasker/canu/masked/01_canu_hybridscaffold_T8.quiver5.pbj3.quiver3.pbj3.quiver3.pilon12.fasta.out.gff
GFF2=/home/jurban/scratch/repeatmasker/canu/alternative/canu_fams/01_canu_hybridscaffold_T8.quiver5.pbj3.quiver3.pbj3.quiver3.pilon12.fasta.out.gff
DB=../../../../protein/db/protein_for_maker_joined_names
PROT=../../../../protein/protein_for_maker_joined_names.fasta

source ~/software/maker/source-maker.sh
blastimap2.1 ${DB} ${REP} -bx -culling_limit 5 -qcov_hsp_perc 1 -evalue 1e-5 > results.paf 
cut -f 6 results.paf | sort | uniq | awk '{gsub(/\\/,"\\\\");gsub(/\|/,"\\|"); gsub(/\(/,"\\("); gsub(/\)/,"\\)"); gsub(/\[/,"\\["); gsub(/\]/,"\\]"); gsub(/\./,"\\."); print}' > regex.txt 
extractFastxEntries.py -R regex.txt -f ${PROT} > results.prot.fasta

filterpaf.py -M -MR 500,500,100,100,1 --same_strand --strandsort -i results.paf > merged.paf

cut -f 1 results.paf | tr "#" "\t" | cut -f 1 | sort | uniq > uniq_rep_fams.txt
cut -f 1 results.paf | sort | uniq > uniq_rep_fams.regex.txt
extractFastxEntries.py -R uniq_rep_fams.regex.txt -f ${REP} > uniq_rep_fams.fasta

awk '$1 !~ /#/ {OFS="\t"; gsub(/Motif:|"/,""); gsub(/\ /,"\t"); print}' ${GFF} | grep.py -p uniq_rep_fams.txt -C 1 -f - -c 10 > uniq_rep_fams.gff
awk '$1 !~ /#/ {OFS="\t"; gsub(/Motif:|"/,""); gsub(/\ /,"\t"); print}' ${GFF2} | grep.py -p uniq_rep_fams.txt -C 1 -f - -c 10 > uniq_rep_fams_2.gff

awk '{OFS="\t"; a[$10]+=1; b[$10]+=$5-$4+1}END{for (e in a) print e,a[e],b[e],b[e]/a[e]}' uniq_rep_fams.gff > uniq_rep_fams.abundance.txt
collapse.py -i results.paf -c 1 -c2 6 -o distinct --strings | tr "#" "\t" > uniq_rep_fams.conversion.txt
grep.py -p uniq_rep_fams.abundance.txt -C 1 -f uniq_rep_fams.conversion.txt -c 1 | paste <( sort -k1,1 uniq_rep_fams.abundance.txt ) <( sort -k1,1  - ) | awk '$1==$5 {OFS="\t"; print $1,$2,$3,$4,$6,$7,$8}' > uniq_rep_fams.abundance-and-annotation.txt

awk '{OFS="\t"; a[$10]+=1; b[$10]+=$5-$4+1}END{for (e in a) print e,a[e],b[e],b[e]/a[e]}' uniq_rep_fams_2.gff > uniq_rep_fams.abundance_2.txt
grep.py -p uniq_rep_fams.abundance_2.txt -C 1 -f uniq_rep_fams.conversion.txt -c 1 | paste <( sort -k1,1 uniq_rep_fams.abundance_2.txt ) <( sort -k1,1  - ) | awk '$1==$5 {OFS="\t"; print $1,$2,$3,$4,$6,$7,$8}' > uniq_rep_fams.abundance-and-annotation_2.txt



Input:
PAF     = results.paf
        = output of blastimap2.2 from blastx of repeatmodeler repeats to protein database.
        = repeat names (column 1) are like "canu_rnd-1_family-11#DNA" and "canu_rnd-1_family-141#Unknown"
        = protein blast database should be made to include entire informative names (connect white space with _ before making)

MERPAF  = merged.paf
        = repeat names (column 1)
        
GFF     = uniq_rep_fams.gff
        = repeat names (column 10 when broken on whitespace)

GFF2     = uniq_rep_fams_2.gff
        = repeat names (column 10 when broken on whitespace)

TXT     = uniq_rep_fams.abundance-and-annotation.txt
        = repeat names (column 1)

TXT     = uniq_rep_fams.abundance-and-annotation_2.txt
        = repeat names (column 1) 

FASTA   = results.prot.fasta
        = protein names

REPFASTA    = uniq_rep_fams.fasta
        

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--PAF',
                   type=str, default='results.paf',
                   help='''PAF''')
parser.add_argument('--MERPAF',
                   type=str, default='merged.paf',
                   help='''PAF''')
parser.add_argument('--GFF',
                   type=str, default='uniq_rep_fams.gff',
                   help='''GFF''')
parser.add_argument('--GFF2',
                   type=str, default='uniq_rep_fams_2.gff',
                   help='''GFF''')
parser.add_argument('--TXT',
                   type=str, default='uniq_rep_fams.abundance-and-annotation.txt',
                   help='''TXT''')
parser.add_argument('--TXT2',
                   type=str, default='uniq_rep_fams.abundance-and-annotation_2.txt',
                   help='''TXT2''')
parser.add_argument('--FASTA',
                   type=str, default='results.prot.fasta',
                   help='''FASTA''')
parser.add_argument('--REPFASTA',
                   type=str, default='uniq_rep_fams.fasta',
                   help='''FASTA''')

parser.add_argument('--OUT',
                   type=str, default='manually-filtered-uniq_rep_fams.abundance-and-annotation',
                   help='''Prefix for output text files. Will make pre.marked_as_repeats.txt and pre.marked_as_genes.txt.
Default: manually-filtered-uniq_rep_fams.abundance-and-annotation''')

parser.add_argument('--yes_to_all',
                   action='store_true', default=False,
                   help='''Use this to blast through without user input. Mainly for testing.''')

parser.add_argument('--prior_hints',
                   type=str, default=False,
                   help='''Directory containing previous input PAF and output txt files. Input paf should be "results.paf".
Output txt files should start with "manually-filtered-uniq_rep_fams.abundance-and-annotation" by default.
Moreover, there should be 2 files with that prefix ending in ".marked_as_repeats.txt" and ".marked_as_genes.txt".''')

parser.add_argument('--repeat_mapping_hints',
                   type=str, default=False,
                   help='''If working with a new repeat library, you can make a mapping between the new and old libraries, then use some information from the old to help decisions with new.
Point this to the repeat to repeat mapping PAF.
The new repeat names should be column 1 (query) and old column 6 (target).
This option requires --prior_hints to also be used.''')


args = parser.parse_args()

if args.repeat_mapping_hints:
    try:
        assert args.prior_hints is not False
    except AssertionError:
        print "Must also use --prior_hints when using --repeat_mapping_hints"
        assert args.prior_hints is not False


## Read in PAF:
paf = defaultdict(list)
with open(args.PAF) as f:
    for line in f:
        line = line.strip()
        name = line.split()[0].split('#')[0]
        paf[name].append(line)
merpaf = defaultdict(list)
with open(args.MERPAF) as f:
    for line in f:
        line = line.strip()
        name = line.split()[0].split('#')[0]
        merpaf[name].append(line)
pafheader = ['repeat','rlen','rstart','rend','strand','hit','hlen','hstart','hend','matches','alnlen','bitscore','NM','PI','EV']

## Read in GFF
gff = defaultdict(list)
with open(args.GFF) as f:
    for line in f:
        line = line.strip()
        name = line.split()[9]
        gff[name].append(line)
gff2 = defaultdict(list)
with open(args.GFF2) as f:
    for line in f:
        line = line.strip()
        name = line.split()[9]
        gff2[name].append(line)
gffheader = ['contig','program','class','start','end','score','strand','.','target','repeat','rstart','rend']

## Read in TXT
txt = defaultdict(list)
with open(args.TXT) as f:
    for line in f:
        line = line.strip()
        name = line.split()[0]
        txt[name].append(line)
txt2 = defaultdict(list)
with open(args.TXT2) as f:
    for line in f:
        line = line.strip()
        name = line.split()[0]
        txt2[name].append(line)
txtheader = ['repeat', 'count', 'span', 'mean-span', 'status', 'npaflines', 'matches']

## Read in FASTA:
fasta = {}
for fa in SeqIO.parse(args.FASTA, "fasta"):
    fasta[str(fa.id)] = str(fa.seq)


## Read in Repeat FASTA
repfasta = {}
for fa in SeqIO.parse(args.REPFASTA, "fasta"):
    repfasta[str(fa.id)] = str(fa.seq)

def get_prot_seq(hit):
    try:
        return fasta[hit]
    except KeyError:
        for key in fasta.keys():
            if hit in key or key in hit:
                return fasta[key]
    return "HIT MISSING FROM FASTA - DEBUG"


def get_rep_seq(hit):
    try:
        return repfasta[hit]
    except:
        pass
    
    try:
        hlen = len(hit)
        seqs=''
        for key in repfasta.keys():
            klen = len(key)
            gate = (klen < hlen and key == hit[:klen]) or (hlen < klen and key[:hlen] == hit)
            if gate:
                seqs+='>'+key+'\n'+repfasta[key]+'\n' 
        if not seqs:
            assert KeyError
        else:
            return seqs
                
    except KeyError: ## this is a band aid
        seqs=''
        for key in repfasta.keys():
            if hit in key or key in hit:
                seqs+='>'+key+'\n'+repfasta[key]+'\n'        
        return seqs
    return "HIT MISSING FROM REPFASTA - DEBUG"

def get_prot_seqs():
    hits = set([])
    for pafline in paf[repeat]:
        hit = pafline.split()[5]
        hits.add(hit)
    for hit in sorted(list(hits)):
        print ">"+hit
        print get_prot_seq(hit)
    print

def get_n_prot_seqs():
    hits = set([])
    for pafline in paf[repeat]:
        hits.add(pafline.split()[5])

    return len(hits)

def get_rep_seqs():
    print ">"+repeat
    print get_rep_seq(repeat)
    print


def get_info():
    #print "#"*20
    npaf = len(paf[repeat])
    npaf2 = len(merpaf[repeat])
    ngff = len(gff[repeat])
    ngff2 = len(gff2[repeat])
    nprot = get_n_prot_seqs()
    
    print repeat
    print npaf, "PAF line(s)"
    print npaf2, "MERPAF line(s)"
    print ngff, "GFF line(s)"
    print ngff2, "GFF2 line(s)"
    print nprot, "proteins"
    print
    print '\t'.join(txtheader)
    for txtline in txt[repeat]:
        print txtline

    print
    print '\t'.join(txtheader)
    for txtline in txt2[repeat]:
        print txtline

    print

def get_txt_alt():
    try:
        t1 = txt[repeat][0].strip().split()
    except:
        t1 = ['-']*7
    #['repeat', 'count', 'span', 'mean-span', 'status', 'npaflines', 'matches']
    try:
        t2 = txt2[repeat][0].strip().split()
    except:
        t2 = ['-']*7
    for i in range(len(txtheader)-1):
        print '\t'.join([txtheader[i], t1[i], t2[i]])
    print
    print txtheader[-1]
    prots = list(set(t1[-1].split(',') + t2[-1].split(',')))
    for prot in prots:
        print prot
    print


def get_paf():
    print pafheader
    for pafline in paf[repeat]:
        print pafline
    print

def get_merpaf():
    print pafheader
    for pafline in merpaf[repeat]:
        print pafline
    print

def get_gff():
    print gffheader
    for gffline in gff[repeat]:
        print gffline
    print

def get_gff2():
    print gffheader
    for gffline in gff2[repeat]:
        print gffline
    print

def gene_it():
    ltxt = len(txt[repeat])
    ltxt2 = len(txt2[repeat])
    if ltxt2 >= 1:
        for txtline in txt2[repeat]:
            mark_as_gene.write(txtline+extra_info+'\n')
    elif ltxt >= 1:
        for txtline in txt2[repeat]:
            mark_as_gene.write(txtline+extra_info+'\n')
    else:
        mark_as_gene.write(repeat+ ' ' + 'DEBUG_NEEDED' + extra_info+'\n')
def repeat_it():
    ltxt = len(txt[repeat])
    ltxt2 = len(txt2[repeat])
    if ltxt2 >= 1:
        for txtline in txt2[repeat]:
            mark_as_repeat.write(txtline+extra_info+'\n')
    elif ltxt >= 1:
        for txtline in txt2[repeat]:
            mark_as_repeat.write(txtline+extra_info+'\n')
    else:
        mark_as_repeat.write(repeat+ ' ' + 'DEBUG_NEEDED' + extra_info+'\n')

## DEBUGGING:
##print '\n'.join(paf['canu_rnd-5_family-807'])
##print '\n'.join(gff['canu_rnd-5_family-807'])
##print '\t'.join(txtheader)
##print '\n'.join(txt['canu_rnd-5_family-807'])
##for pafline in paf['canu_rnd-5_family-807']:
##    hit = pafline.split()[5]
##    print hit, get_seq(hit)

def get_raw_input(msg):
    ans = raw_input(msg)
    while ans not in ('y','yes','n','no'):
        print ans, 'is not a viable response....'
        print 'try:','y','yes','n','no'
        ans = raw_input(msg)
    return ans

def clear():
    print
    print '#'*100
    print





def get_repeats_already_processed():
    rephints = defaultdict(list)
    alreadydone = set([])
    print "....already processed....."
    cat = ['REPEAT','GENE']
    catd = defaultdict(int)
    i=-1
    for fh in (args.OUT+'.marked_as_repeats.txt', args.OUT+'.marked_as_genes.txt'):
        i+=1
        with open(fh) as f:
            for line in f:
                line = line.strip().split('\t')
                print line[0], cat[i], line[-1]
                catd[cat[i]] += 1
                alreadydone.add(line[0])
                rephints[line[0]] = [cat[i], line[-1]]
    N = len( list(set(txt.keys() + txt2.keys())) )
    DONE = len(alreadydone)
    print
    print "Already processed", DONE, "repeats of", N
    print catd['REPEAT'], 'repeats'
    print catd['GENE'], 'genes'
    print
    print "There are", N-DONE, "repeats left to do..."
    print 
    ans = 'n'
    while ans in ('n','no'):
        ans = get_raw_input('Shall we begin? y/n. ::: ')
    return alreadydone, rephints

def update_hints(hints, rephints):
    hits = set([])
    for pafline in paf[repeat]:
        hit = pafline.split()[5]
        hits.add(hit)
    for hit in list(hits):
        hints[hit].append( [repeat] + rephints[repeat] )
    return hints, rephints

def get_hints():
    hits = set([])
    for pafline in paf[repeat]:
        hit = pafline.split()[5]
        hits.add(hit)
    totalentries = 0
    totalhits = 0
    totalrepeats = 0
    hits_assoc_with_rep = 0
    hitreps = defaultdict(int)
    hitents = defaultdict(int)
    for hit in list(hits):
        totalhits += 1
        repcnt = 0
        entcnt = 0
        try:
            for entry in hints[hit]:
                totalentries += 1
                entcnt += 1
                out = '\t'.join([hit]+entry)
                print out
                if 'REPEAT' in out:
                    repcnt += 1
                    totalrepeats += 1
            if repcnt > 0:
                hits_assoc_with_rep += 1
            hitreps[hit] = repcnt
            hitents[hit] = entcnt
        except:
            continue
    print
    print totalrepeats, 'of', totalentries, 'previous entries/hints assoc with repeats.'
    print hits_assoc_with_rep, 'of', totalhits, 'current hits assoc with previous repeat hints.'
    print
    print totalrepeats, 'of', totalentries, 'previous entries/hints assoc with repeats.', hits_assoc_with_rep, 'of', totalhits, 'current hits assoc with previous repeat hints.'
    print
    for key in hitents.keys():
        print key, hitreps[key], 'of', hitents[key], 'entries assoc with previous repeat hints.'
    print


def build_prior_hints(pre='manually-filtered-uniq_rep_fams.abundance-and-annotation'):
    # msg
    print "....collecting info on prior hints....."
    # init
    rephints = defaultdict(list)
    rep_assoc_hits = defaultdict(set)
    prior_hints = defaultdict(list)
    cat = ['REPEAT','GENE']
    hits = set([])
    i=-1
    # iter
    for fh in (args.prior_hints+pre+'.marked_as_repeats.txt', args.prior_hints+pre+'.marked_as_genes.txt'):
        i+=1
        with open(fh) as f:
            for line in f:
                line = line.strip().split('\t')
                # Get priorrepeat and resulting comment. NOTE: prior repeats may share names with curr, but have nothing to do wth each other. Careful.
                rephints[line[0]] = [cat[i], line[-1]]
    N = len(rephints.keys())
    print "Found ", N, "rep hints...."
    print
    ## Using prior input PAF, get each hit correspondng to a repeat
    with open(args.prior_hints+'results.paf') as f:
        for line in f:
            line = line.strip()
            repeat = line.split()[0].split('#')[0]
            hit = line.split()[5]
            rep_assoc_hits[repeat].add(hit)
    ## Correlate the hits assoc with a repeat with the repeat result
    for repeat in rephints.keys():
        for hit in list( rep_assoc_hits[repeat] ):
            prior_hints[hit].append(rephints[repeat])
    N = len(prior_hints.keys())
    print "Found ", N, "prior hints on hits in PAF...."
    ans = 'n'
    while ans in ('n','no'):
        ans = get_raw_input('Shall we move on? y/n. ::: ')
    print
    return prior_hints, rephints


def get_prior_hints():
    hits = set([])
    for pafline in paf[repeat]:
        hit = pafline.split()[5]
        hits.add(hit)
    totalentries = 0
    totalhits = 0
    totalrepeats = 0
    hits_assoc_with_rep = 0
    hitreps = defaultdict(int)
    hitents = defaultdict(int)
    for hit in list(hits):
        totalhits += 1
        repcnt = 0
        entcnt = 0
        try:
            for entry in prior_hints[hit]:
                totalentries += 1
                entcnt += 1
                out = '\t'.join([hit]+entry)
                print 'PRIOR INFORMATION:\n', out, '\n'
                if 'REPEAT' in out:
                    repcnt += 1
                    totalrepeats += 1
            if repcnt > 0:
                hits_assoc_with_rep += 1
            hitreps[hit] = repcnt
            hitents[hit] = entcnt
        except:
            continue
    print
    print 'PRIOR INFORMATION:', totalrepeats, 'of', totalentries, 'previous entries/hints assoc with repeats.'
    print 'PRIOR INFORMATION:', hits_assoc_with_rep, 'of', totalhits, 'current hits assoc with previous repeat hints.'
    print
    print 'PRIOR INFORMATION:', totalrepeats, 'of', totalentries, 'previous entries/hints assoc with repeats.', hits_assoc_with_rep, 'of', totalhits, 'current hits assoc with previous repeat hints.'
    print
    for key in hitents.keys():
        print 'PRIOR INFORMATION:', key, hitreps[key], 'of', hitents[key], 'entries assoc with previous repeat hints.'
    print
    
def build_repeat_mapping_hints():
    print 'Building map of current repeats to prior set....'
    maphints = defaultdict(list)
    with open(args.repeat_mapping_hints) as f:
        for line in f:
            line = line.strip()
            name = line.split()[0].split('#')[0]
            maphints[name].append(line)
    print 'Found', len(maphints.keys()), 'current repeats with mappings to older ones.'
    ans = 'n'
    while ans in ('n','no'):
        ans = get_raw_input('Shall we move on? y/n. ::: ')
    print
    return maphints

def get_map_hints():
    for mapping in maphints[repeat]:
        mapping = mapping.strip().split()
        old, oldclass = mapping[5].split('#')
        print 'MAPPING HINT:'
        print '\t'.join(mapping)
        for hint in prior_rephints[old]:
            print hint
        print


def collect(extra_info):
    if extra_info == '\t':
        pre = ''
    else:
        pre = ' '
    if not args.yes_to_all:
        ans = get_raw_input('Would you like to add information to output? y/n  :::  ')
    if not args.yes_to_all and ans.lower() in ('y', 'yes'):
        extra_info += pre + raw_input('Please type in info to be appended to output. Then press enter. \n\n Text here :: ')
    print
    return extra_info

## EXECUTE:
prior_hints = False
maphints = False
if args.prior_hints is not False:
    prior_hints, prior_rephints = build_prior_hints()
if args.repeat_mapping_hints is not False:
    maphints = build_repeat_mapping_hints()
    
hints = defaultdict(list)
rephints = defaultdict(list)
going_msgs = ("You're a trooper!", "Gettin' it done.", "woot!", "Not one to procrastinate, are ya?", "Okay okay okay... I'll pull up the next one... jeeeeez.", "Congrats - you earned my respect.", "Looks like we are gonna need some coffee.", "Let me go grab a chair... this is gonna be a while.","Excellent! Good for you. (You know you can pick back up where you leave off right?). Terrific.", "Choo choooooooo!", "You're not gonna let this go, are you?", "You and me baby, we ain't nothin' but mammals, so let's do it like they do on the discovery channel.", "You said yes. Again.", "Here we go again.", "From all of us here inside your computer: we are unionizing.", "...it's been a long day's night and i've been working like a dog...", "Work it.", "Yeah!", "I thought you'd say that.", "Here's your prize:", "And the next one is:", "Uh huh.")
len_going_msgs = len(going_msgs)
ans = 'n'
if os.path.isfile(args.OUT+'.marked_as_repeats.txt') and os.path.isfile(args.OUT+'.marked_as_genes.txt'):
    ans = get_raw_input('Would you like to pick back up where you left off? y/n  :::  ')
if not args.yes_to_all and ans.lower() in ('y', 'yes'):
    continuing = True
    alreadydone, rephints = get_repeats_already_processed()
    mark_as_repeat = open(args.OUT+'.marked_as_repeats.txt', 'a')
    mark_as_gene = open(args.OUT+'.marked_as_genes.txt', 'a')
else:
    continuing = False
    alreadydone = set([])
    mark_as_repeat = open(args.OUT+'.marked_as_repeats.txt', 'w')
    mark_as_gene = open(args.OUT+'.marked_as_genes.txt', 'w')

ans = 'yes'
clear()
NREPS = len( list(set(txt.keys() + txt2.keys())) )
try:
    for repeat in sorted( list(set(txt.keys() + txt2.keys())) ):
        extra_info = '\t'
        
        if continuing:
            if repeat in alreadydone:
                hints, rephints = update_hints(hints, rephints)
                continue
          
        get_info()

        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see txt/txt2 summary in alt format? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_txt_alt()
            extra_info = collect(extra_info)

        if not args.yes_to_all and maphints:
            ans = get_raw_input('Would you like to see hints from corresponding repeats in PRIOR library (if any exist)? y/n  :::  ')
        if maphints and ans.lower() in ('y', 'yes'):
            get_map_hints()
            extra_info = collect(extra_info)

        if not args.yes_to_all and prior_hints:
            ans = get_raw_input('Would you like to see PRIOR hints if any exist? y/n  :::  ')
        if prior_hints and ans.lower() in ('y', 'yes'):
            get_prior_hints()
            extra_info = collect(extra_info)
        
        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see hints if any exist? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_hints()
            extra_info = collect(extra_info)
            
        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see PAF lines? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_paf()

        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see MERGED PAF lines? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_merpaf()

        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see GFF lines? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_gff()

        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see GFF lines in GFF2? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_gff2()

        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see repeat seq? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_rep_seqs()

        if not args.yes_to_all:
            ans = get_raw_input('Would you like to see protein seqs? y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
            get_prot_seqs()

##        if not args.yes_to_all:
##            ans = get_raw_input('Would you like to add information to output? y/n  :::  ')
##        extra_info = ''
##        if not args.yes_to_all and ans.lower() in ('y', 'yes'):
##            extra_info += '\t' + raw_input('Please type in info to be appended to output. Then press enter. \n\n Text here :: ')
##        print
        extra_info = collect(extra_info)
        
        if not args.yes_to_all:
            ans = get_raw_input('Would you call this a gene? \n i.e. should this be left un-masked during annotation? \n\t\t y/n  :::  ')
        if ans.lower() in ('y', 'yes'):
             gene_it()
             rephints[repeat] = ['GENE', extra_info]
        else:
            repeat_it()
            rephints[repeat] = ['REPEAT', extra_info]
        hints, rephints = update_hints(hints, rephints)

        clear()

        NDONE = len(rephints.keys())
        print "Processed", NDONE, "repeats of", NREPS, '...'
        print "There are", NREPS-NDONE, "repeats left to process..."
        
        if not args.yes_to_all:
            ans = get_raw_input('Would you like to keep going? y/n  :::  ')
        extra_info = ''
        if ans.lower() in ('n', 'no'):
            print "You can pick back up where you left off later. Good bye! :)"
            break
        else:
            print choice(going_msgs)+"\n"
            sleep(1)
except KeyboardInterrupt:
    print
    print
    print "I see you impatiently canceled... Well I will save progress anyway...."
    print
    sleep(1)

mark_as_gene.close()      
mark_as_repeat.close()



if not os.path.exists('manual_repeat_filtering_backup'):
    os.mkdir('manual_repeat_filtering_backup')
now = '-'.join(''.join('_'.join(str(datetime.datetime.now()).split()).split('.')[:-1]).split(':'))
for e in ('repeats','genes'):
    print
    cmd = 'cp ' + args.OUT+'.marked_as_' + e + '.txt manual_repeat_filtering_backup/' + args.OUT+'.marked_as_'+ e + '.' + now + '.txt' 
    print cmd
    os.system(cmd)
print
