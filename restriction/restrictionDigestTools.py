import sys, re, string
import argparse
import random
import os
from itertools import product
from Bio import SeqIO
from collections import defaultdict
from restrictionEnzymes import *

####################################################################################################
'''---------------------------[ Debug functions ]----------------------------'''
####################################################################################################

def debug(msg, args):
    if args.debug or args.verbose:
        sys.stderr.write('[::restrictionDigest::] ' + msg +'\n')






##def fastaToPython(fastaFile):
##    seqs = open(fastaFile, 'r')
##    fasta = seqs.readlines()
##    seqs.close()
##    headers = []
##    DNAstrings = []
##    index=0
##    ## puts names in 'headers'
##    ## puts sequence in 'DNAstrings'
##    ##      with same index as corresponding name in 'headers'
##    ## To do this, I combine all lines in between '>' header lines
##    ## Moreover, I make into a single string not separated by '\n'
##    for element in fasta:
##        if element[0] == '>':
##            if element[-1] == '\n':
##                headers = headers + [element[:-1]]
##            else:
##                headers = headers + [element]
##            DNAstrings = DNAstrings + ['']
##            index+=1
##        else:
##            if element[-1] == '\n':
##                DNAstrings[index-1] = DNAstrings[index-1] + element[:-1]
##            else:
##                DNAstrings[index-1] = DNAstrings[index-1] + element
##    
##    return headers, DNAstrings


##def fastaToPython(fastaFile, seqdict={}):
##    with open(fastaFile, 'r') as seqs:
##        fasta = seqs.readlines()
##    ## puts names in 'headers'
##    ## puts sequence in 'DNAstrings'
##    ##      with same index as corresponding name in 'headers'
##    ## To do this, I combine all lines in between '>' header lines
##    ## Moreover, I make into a single string not separated by '\n'
##    for element in fasta:
##        if element[0] == '>':
##            header = element.strip()
##            debug("Loading : " + header )
##            seqdict[header] = ''
##        else:
##            seqdict[header] += element.strip()
##    return seqdict

def fastaToPython(fastaFile, seqdict={}):
    with open(fastaFile, 'r') as fasta:
        for element in fasta:
            if element[0] == '>':
                header = element.strip()
                debug("Loading : " + header )
                seqdict[header] = ''
            else:
                seqdict[header] += element.strip()
    return seqdict

##def allACGT(x):
##    xlen = len(x)
##    n = sum(e in 'ACGT' for e in x)
##    return n == xlen
##

def rseqify(rseq):
    ''' given seq return all possible seqs'''
    d = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'N':'ACGT'}
    return [''.join(j) for j in product(*(d[b] for b in rseq) ) ] 


def regexify(rseq):
    ''' given seq return all possible seqs'''
    d = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'N':'ACGT'}
    regex = ''
    oldb = rseq[0]
    cnt = 1
    for i in range(1,len(rseq)):
        b = rseq[i]
        if b == oldb:
            cnt += 1
            continue

        
        if len(d[oldb]) == 1:
            regex += oldb
        else:
            regex += '[' + d[oldb] + ']'

        if cnt > 1:
            regex += '{' + str(cnt) + '}'
        cnt = 1
        oldb = b

    ## PROCESS LAST ONE
    if len(d[b]) == 1:
        regex += b
    else:
        regex += '[' + d[b] + ']'

    if cnt > 1:
        regex += '{' + str(cnt) + '}'

    return regex


def revcomp(rseq):
    intab=  'actguACTGUNRY'
    outtab= 'tgacaTGACANYR'
    transtab = string.maketrans(intab, outtab)
    return rseq.translate(transtab)[-1::-1]    

##def regexifyrevcomp(rseq):
##    rc = revcomp(rseq)
##    if rc == rseq:
##        return regexify(rseq)
##    return regexify(rc)

def rseqify_list(x):
    ''' given list of restriction seqs, decode if needed for N,R,Y'''
    rseqs = []
    for rseq in x:
        rseqs += rseqify(rseq)
    return rseqs
    

def rseqify_list2(x):
    ''' given list of restriction enzymes, decode if needed for N,R,Y'''
    enames = []
    rseqs = []
    for e in x:
        seqs = rseqify(enzymes[e])
        lenseqs = len(seqs)
        if lenseqs == 1:
            enames.append( e )
        elif lenseqs > 1:
            enames += [e+"_"+str(i) for i in range(lenseqs)]
        rseqs += seqs

        #RCs
        rc = revcomp(enzymes[e])
        if rc != enzymes[e]:
            rcseqs = rseqify(rc)
            lenrcseqs = len(rcseqs)
            if lenrcseqs == 1:
                enames.append( e + '_rc' )
            elif lenrcseqs > 1:
                enames += [e+"_rc_"+str(i) for i in range(lenrcseqs)]
            rseqs += rcseqs
                
    rsitelengths = [len(e) for e in rseqs]
    return enames, rseqs, rsitelengths  


def regexify_list(x):
    ''' given list of restriction enzymes, decode if needed for N,R,Y'''
    enames = []
    rseqs = []
    rsitelengths = []
    for e in x:
        enames.append( e )
        rseqs.append( regexify(enzymes[e]) )
        rsitelengths.append( len(enzymes[e]) )

        #RCs
        rc = revcomp(enzymes[e])
        if rc != enzymes[e]:
            enames.append( e + '_rc' )
            rseqs.append( regexify(rc) )
            rsitelengths.append( len(rc) )
    return enames, rseqs, rsitelengths


## THIS WAS RETURNING 1-based starts....
##def substringPositions1(DNAstring, DNAsubstring):
##    """returns 1-based positions"""
##    window = len(DNAsubstring)
##    locations = []
##    for i in range(len(DNAstring)-window+1):
##        if DNAstring[i:i+window] == DNAsubstring:
##            locations = locations + [i+1] ## only diff from 0-based fxn
##    return locations


## NOW RETURNS 0-BASED STARTS
def substringPositions1(DNAstring, DNAsubstring):
    """returns 0-based positions"""
    window = len(DNAsubstring)
    locations = []
    for i in range(len(DNAstring)-window+1):
        if DNAstring[i:i+window] == DNAsubstring:
            locations = locations + [i] ## only diff from 0-based fxn
    return locations


def regexStartPositions(DNAstring, DNAsubstring):
    ''' Has same outout as 0-based substringPositions1 fxn'''
    regex = re.compile(DNAsubstring)
    locations = []
    for m in re.finditer(regex, DNAstring):
        locations.append(m.start())
    return locations





def seqlist2regexlist(l):
    ''' given list of sequences, return list of regexes'''
    return 

def regexlist2regex(l):
    ''' given list of sequences, return list of regexes'''
    return '|'.join(l)

def bedify_startslist(header, starts, motif):
    bedlist = []
    motiflen = len(motif)
    for start in starts:
        bedlist.append( (header, start, start+motiflen, motif) )
    return bedlist

def calculateLengths(starts, seq, lengths=[], beds=[], getbed=False, addXtoBedEnd=[], labels=[]):
    num = len(starts)
    #debug("calculateLengths Fxn: "+str(num)+" "+str(len(lengths))+" "+str(len(beds)))
    if num == 1:
        #debug("calculateLengths Fxn: 1")
        lengths += [starts[0]]
        lengths += [len(seq)-starts[0]] 
        if getbed:
            beds.append( (0, starts[0]+addXtoBedEnd[0], '5p:'+labels[0]) )
            beds.append( (starts[0], len(seq), labels[0]+':3p' ) )
    if num > 1:
        #debug("calculateLengths Fxn: MORE THAN 1")
        lengths += [starts[0]]
        beds.append( (0, starts[0]+addXtoBedEnd[0], '5p:'+labels[0]) )
        for i in range(1, len(starts)):
            lengths += [starts[i] - starts[i-1]] ## Length is defined from start to start 
            if getbed:
                beds.append( (starts[i-1], starts[i]+addXtoBedEnd[i], labels[i-1]+':'+labels[i]) )
        if len(starts) != starts[i-1]:
            lengths += [len(seq)-starts[i]] 
            if getbed:
                beds.append( (starts[i], len(seq), labels[i]+':3p') )
    return lengths, beds



## Beta1 version
##def calculateCircularLengths(starts, seq, lengths=[], beds=[]):
##    """Need to be able to process plasmids"""
##    num = len(starts)
##    if num == 1:
##        lengths += [len(seq)]
##    if num > 1:
##        new1spot = starts[0]
##        seq = seq[new1spot :] + seq[:new1spot]
##        for i in range(len(starts)):
##            starts[i] = starts[i] - new1spot + 1
##        for i in range(1, len(starts)):
##            lengths += [starts[i] - starts[i-1]]
##        if len(starts) != starts[i-1]:
##            lengths += [len(seq)-starts[i]+1] ## I think +1 is wrong
##    return lengths


def calculateCircularLengths(starts, seq, lengths=[], beds=[], getbed=False, addXtoBedEnd=[], labels=[]):
    """Need to be able to process plasmids"""
    num = len(starts)
    if num == 1:
        lengths += [len(seq)]
        if getbed:
            beds.append( (starts[0], len(seq)+starts[0]+addXtoBedEnd[0], labels[0]+':'+labels[0]) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
    if num > 1:
        for i in range(1, len(starts)): # skip first start: starts[0]
            lengths += [starts[i] - starts[i-1]]
            if getbed:
                beds.append( (starts[i-1], starts[i]+addXtoBedEnd[i], labels[i-1]+':'+labels[i]) )
        # Last one = len-to-end + 0-to-start[0]
        lengths += [len(seq)-starts[i] + starts[0]] 
        if getbed:
            beds.append( (starts[i], len(seq)+starts[0]+addXtoBedEnd[i] , labels[i]+':'+labels[0]) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
    return lengths, beds

##def calculateCircularLengths(starts, seq, lengths=[], beds=[]):
##    """Need to be able to process plasmids"""
##    num = len(starts)
##    if num == 1:
##        lengths += [len(seq)]
##        beds.append( (starts[0], len(seq)+starts[0]) ) ## difficult situation.... for BED on circle... I decided to double the circle here.
##    if num > 1:
##        new1spot = starts[0]
##        seq = seq[new1spot :] + seq[:new1spot]
##        for i in range(len(starts)):
##            starts[i] = starts[i] - new1spot + 1
##        for i in range(1, len(starts)):
##            lengths += [starts[i] - starts[i-1]]
##            beds.append( (starts[i-1], starts[i]) )
##        if len(starts) != starts[i-1]:
##            lengths += [len(seq)-starts[i]+1]
##    return lengths

def identity(x):
    return x


class Digest(object):
    def __init__(self, header, seq, enzymenames, recseqs):
        '''
        header is single fasta header or other string
        seq is single associated fasta sequence string
        names is a list of N strings for enzyme names or naming of recseqs somehow
        recseqs is a list of N associated strings/recognition sequences (can do single enzymedict in future, but code is made for comma-sep list atm)
        
        '''
        self.header = seqname
        self.seq = seq
        self.names = enzymenames
        self.recseqs = enzymeseqs
        self.n = len(enzymenames)
        assert self.n == len(enzymeseqs)
        self.niter = range(self.n)
        # Get Start positions of each enzymeseq
        self.startsList = [substringPositions1(self.seq, self.recseqs[i]) for i in self.niter]
        




def check():
    nseq=len(seqs)
    nhead = len(headers)
    debug('Loaded ' + str(nhead) + ' headers and ' + str(nseq) + ' associated sequences.')
    assert nseq == nhead

    







##################################################################################################
#################################################################################################
##################################################################################################
#################################################################################################
##################################################################################################
#################################################################################################
##################################################################################################
#################################################################################################
##################################################################################################
#################################################################################################
##################################################################################################
#################################################################################################
##################################################################################################
#################################################################################################




            
def VER2PIPELINE(header, seq, args, rseq_translation_fxn, start_positions_fxn):
    ##################################################################################################
    '''---------------------------[ Early processing of arguments ]----------------------------'''
    ##################################################################################################
    debug('Stage: Early processing of arguments.',args)
    debug('Setting linear vs. circle.',args)
    linear=True
    if args.circular:
        linear = False
    ##########################################################################
    '''---------------------------[ Execution ]----------------------------'''
    ##########################################################################
    debug('Stage: Execution.',args)



    ##############################################################################################
    '''---------------------------[ PROCESS RESTRICTION SITES ]----------------------------'''
    ##############################################################################################
    debug('Stage: Processing restriction sites.',args)

    ## Process restriction site(s) as list
    if args.sites_processed:
        pass
    else:
        args.sites_processed = True
        if args.restrictionEnzyme:
            args.restrictionEnzyme, args.restrictionSite, args.restrictionSiteLengths = rseq_translation_fxn([e for e in args.restrictionEnzyme.split(',')])
            #enzymes = {args.restrictionEnzyme[i]:args.restrictionSite[i] for i in range(len(args.restrictionSite))}
            for i in range(len(args.restrictionSite)):
                enzymes[args.restrictionEnzyme[i]] = args.restrictionSite[i]
                enzymelookup[args.restrictionSite[i]].append(args.restrictionEnzyme[i])
            if len(args.restrictionEnzyme) <= 10:
                debug(str(zip(args.restrictionEnzyme, args.restrictionSite)),args)
        else:
            args.restrictionSite = args.restrictionSite.split(',') 
            args.restrictionEnzyme = [e for e in args.restrictionSite]
            for e in args.restrictionEnzyme:
                enzymes[e]=e
                enzymelookup[e].append(e)

        
    ##############################################################################################
    '''---------------------------[ FIND RESTRICTION SITES ]----------------------------'''
    ##############################################################################################
    debug('Stage: Finding restriction sites.',args)

    debug('Identifying restriction site starts in sequence: ' + str(header) ,args)
    if len(args.restrictionSite) > 1:
        debug('More than one enzyme detected - i.e. user provided more than one enzyme or recognition sequence.',args)
        
        # startsList is a list of lists
        startsList = [start_positions_fxn(seq, args.restrictionSite[i]) for i in range(len(args.restrictionSite))]

        n_startsList = len(startsList)
        n_ns_total = 0
        for i in range(n_startsList):
            n_sl_i = len(startsList[i])
            n_ns_total += n_sl_i
            debug('Detected ' + str(n_sl_i) + ' starts for ....' + args.restrictionSite[i],args)
        debug('Detected ' + str(n_ns_total) + ' sites overall.',args)
        debug('Non-empty starts lists came from:',args)
        if args.debug:
            for i in range(len(startsList)):
                if len(startsList[i])>0:
                    debug('\t'+args.restrictionSite[i],args)
                
        #bedlist = list of lists of tuples
        if args.bed:
            pass
            #bedlist = [(header, startsList[i][j], startsList[i][j]+len(args.restrictionSite[i]), args.restrictionSite[i]) for j in range(len(startsList)) for i in range(len(args.restrictionSite))]

        if args.digestTogether:
            if n_ns_total > 0: #Then combine list of lists into single list of starts and sort (as if only a single restrictionSite used)
                debug('DNA will be digested by more than one enzyme together in the same tube.',args)
                #starts = sorted( [x for e in startsList for x in e] )
                #starts = sorted( [x for i in range(len(startsList)) for x in startsList[i] ] )
                ### starts, starts_labels, sitelengths = zip(*sorted( [(x,args.restrictionEnzyme[i], len(args.restrictionSite[i])) for i in range(n_startsList) for x in startsList[i] ], key = lambda x: x[0] )) ## * means to unzip
                starts, starts_labels, sitelengths = zip(*sorted( [(x,args.restrictionEnzyme[i], args.restrictionSiteLengths[i]) for i in range(n_startsList) for x in startsList[i] ], key = lambda x: x[0] )) ## * means to unzip
                debug('Repeat: Detected ' + str(len(starts)) + ' sites overall.',args)
            else:
                starts, starts_labels, sitelengths = [], [], []


    elif len(args.restrictionSite) == 1:
        debug('Only one enzyme detected - i.e. user provided a single enzyme or recognition sequence.',args)

        # if only one restrictionSite, then just use a list of Starts
        starts = start_positions_fxn(seq, args.restrictionSite[0])
        starts_labels = [args.restrictionEnzyme[0]]*len(starts)
        sitelengths = [args.restrictionSiteLengths[0]]*len(starts)
        debug('Detected ' + str(len(starts)) + ' sites overall.',args)




    ##############################################################################################
    '''---------------------------[ CALCULATE LENGTHS ]----------------------------'''
    ##############################################################################################
    debug('Stage: Calculate lengths.',args)


    if args.digestTogether or len(args.restrictionSite) == 1: ## starts is a list
        debug('Single tube reaction.',args)
        if linear:
            debug('Calculating lengths given a linear molecule.',args)
            lengths, intervals = calculateLengths(starts, seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)
        elif args.circular:
            debug('Calculating lengths given a circular molecule.',args)
            lengths, intervals = calculateCircularLengths(starts, seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)

         
    else: ## starts is a list of lists
        ## Get lengths list for each list in list of lists
        ## Store as list of lists
        debug('More than one reaction.',args)
        lengthsList = []
        intervals = []
        for i in range(len(startsList)):
            starts = startsList[i]
            starts_labels = [args.restrictionEnzyme[i]]*len(starts)
            sitelengths = [len(args.restrictionSite[i])]*len(starts)
            if linear:
                debug('Calculating lengths given a linear molecule for reaction ' + str(i) + '.',args)
                lengths, beds = calculateLengths(startsList[i], seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)
                lengthsList += [lengths]
                intervals += [beds]
            elif args.circular:
                debug('Calculating lengths given a circular molecule for reaction ' + str(i) + '.',args)
                lengths, beds = calculateCircularLengths(startsList[i], seq, lengths=[], beds=[], getbed=args.bed, addXtoBedEnd=sitelengths, labels=starts_labels)
                lengthsList += [lengths]
                intervals += [beds]

        if args.loadTogether:
            debug('Loading restriction reactions into same well for gel electrophoresis.',args)
            ## then lengths is list of lengths you would see in a single lane
            ## combine the list of lists into a single list of lengths

            lengths = [x for e in lengthsList for x in e]
            #lengths = sorted( [x for e in startsList for x in e] )
        else:
            debug('NOT loading restriction reactions into same well for gel electrophoresis.',args)
            lengths = []


    ## ADD LENGTH OF UNDIGESTED MOLECULE?
    if args.length and linear:
        debug('Intact linear molecule will be present.',args)
        lengths += [len(seq)]
    elif args.length and args.circular:
        debug('Intact circular molecule(s) will be present.',args)
        lengths += ["super-coiled and nicked circles"]
        

    ## SORT LENGTHS INTO NEW LIST
    debug('Sorting lengths.',args)
    sortedLengths = sorted(lengths, reverse=True)



    ##############################################################################################
    '''---------------------------[ OUTPUTS ]----------------------------'''
    ##############################################################################################
    debug('Stage: Output.',args)

    ## For when >1 site given, but not to be put together
    ## Calculate separate reactions...
    if args.bed:
        if len(args.restrictionSite) == 1 or args.digestTogether:
            for i in range(len(intervals)):
                interval = intervals[i]
                label = interval[2]
                labelseqs = ':'.join([enzymes[e] for e in label.split(':')])
                start_to_start_length = lengths[i]
                interval_length = interval[1]-interval[0]
                bed = [header, interval[0], interval[1], label+'_'+labelseqs, interval_length, '+', start_to_start_length]
                print ('\t').join([str(e) for e in bed])
        else:
            debug('multiple-sites for BED still in dev. Need different solutions for digest together and sep... if doing sep, then just analyze one at a time... duh',args)
            quit()

    elif len(args.restrictionSite) > 1 and not (args.loadTogether or args.digestTogether):
        debug('Lengths are in single well of gel.',args)
        for i in range(len(args.restrictionSite)):
            print args.restrictionSite[i]
            if args.proportion:
                debug('Calculating proportions for each fragment.',args)

                if args.length: ## This may seem redundant with section: ADD LENGTH OF UNDIGESTED MOLECULE? -- but its not b/c that is on "lengths" and this is using "lengthsList below
                    print args.mass*len(seq)/float(len(seq))

                for j in range(len(lengthsList[i])):
                    print args.mass*lengthsList[i][j]/float(len(seq))

                if i < len(args.restrictionSite)-1: 
                    print ## Print space between current and next output.
            else:
                debug('Showing lengths of each fragment as seen in gel.',args)

                if args.length: ## This may seem redundant with section: ADD LENGTH OF UNDIGESTED MOLECULE? -- but its not b/c that is on "lengths" and this is using "lengthsList below
                    print len(seq) 
                    
                for j in range(len(lengthsList[i])):
                    print lengthsList[i][j]

                if i < len(args.restrictionSite)-1: 
                    print ## Print space between current and next output.
    else:
        debug('Lengths are in separate wells of gel for different reactions.',args)
        if args.ordered:
            L = lengths
        else:
            L = sortedLengths
            
        if args.rmap:
            print '\t'.join([str(e) for e in [header, len(seq), len(L)] + L])
        else:
            print ">"+header
            for i in range(len(L)):
                if args.proportion:
                    print args.mass*L[i]/float(len(seq))
                else:
                    print L[i]
    debug('Done with: '+header,args)
    debug('',args)





def ver2pipeline(header, seq, args):
    VER2PIPELINE(header,seq, args, rseqify_list2, substringPositions1)


## TODO: I am not sure if I'd prefer the regex as the BED label or the actual seq
def regex_ver2pipeline(header, seq, args):
    ## IN THIS PIPELINE, EACH ENZYME WILL USED IN AN INDEPENDENT REGEX
    ##  IN regex_ver3pipeline A SINGLE REGEX WILL BE USED AND SITE LABELS WILL BE DEFINED BY THE ENZYMELOOKUP TABLE
    ##   ...A TODO
    VER2PIPELINE(header,seq, args, regexify_list, regexStartPositions)
