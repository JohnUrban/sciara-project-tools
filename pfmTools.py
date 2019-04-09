import random, sys

def addToPfm(seqLen, sequence, posFreqMatDict, currentHeader, seqCount):
    if len(sequence) == seqLen:
        sequence = sequence.upper() ## make all bases same case
        ## for each position along sequence
        for i in range(seqLen):
            ## Take the base at position i
            b_i = sequence[i]
            try:
                ## Add 1 count at position i for b_i = b
                posFreqMatDict[i][b_i] += 1
            except KeyError:
                try: ## perhaps position i is there, but not yet this base
                    posFreqMatDict[i][b_i] = 1
                except KeyError: ## this means position i was not yet there
                    posFreqMatDict[i] = {}
                    posFreqMatDict[i][b_i] = 1
            
    else:
        print "Kill Error: Reached a sequence of different length than previous sequences."
        print "Sequence length established as:", seqLen
        print "Sequence name:", currentHeader
        print "Sequence number:", seqCount
        print "Sequence length:", len(sequence)
        quit()
    return posFreqMatDict




def pfmDictFromFasta(fastaFile):
    """Fasta File should contain N sequences all of the same length.
    If a sequence of a length != to the length of the first sequence analyzed is encountered,
    it will raise an error.
    Note, there should be no empty lines in the fasta file except for the end of file.
    Use: "grep -v ^$ fasta file > new.fa" to fix if you need."""

    ## open connection to the fasta file
    fasta = open(fastaFile, 'r')

    ## initialize a variable to establish sequence length
    seqLen = None

    ## initialize a variable to store current sequences
    sequence = ''

    ## initalize the PFM dictionary
    posFreqMatDict = {}

    ## Do line#1 outside of for loop
    currentHeader = fasta.readline().strip()[1:]
    seqCount = 1

    ## For each line in fasta, determine if it is an empty line (end of file), a header line, or a sequence line and do appropriate actions.
    for line in fasta:
        if line.startswith('>') == False: ## sequence line
            sequence += line.strip() ## Add all but \n at end
        else: ##line starts with '>'
            ## Then there must be a sequence that has just been put together
            ## analyze previous sequence before moving on
            ## 1. Make sure stored sequence is same size as all previous sequences (or establish the seqLen now)
            if seqLen == None:
                ## This should only occur once
                ## establish seqLen 
                seqLen = len(sequence)
                ## establish what the minimum dictionary should be
                for i in range(seqLen):
                    posFreqMatDict[i] = {}
                    for b in 'ACGTN':
                        posFreqMatDict[i][b] = 0
            posFreqMatDict = addToPfm(seqLen, sequence, posFreqMatDict, currentHeader, seqCount)
            currentHeader = line[1:-1]
            seqCount += 1
            sequence = ''
    posFreqMatDict = addToPfm(seqLen, sequence, posFreqMatDict, currentHeader, seqCount)
    fasta.close()
    return posFreqMatDict



def pfmFromPfmDict(posFreqMatDict):
    """Takes in PFM dictionary object -- e.g. output of pfmFromFasta
    Writes it out to file ...."""
    
    ##1. Find the full set of symbols
    bases = set()
    for i in posFreqMatDict:
        for base in posFreqMatDict[i].keys():
            bases.add(base) ## "bases" is a set variable -- thus, if the element 'base' is in the set, it will not be repeated

    ##2. Ensure that each position in dictionary has all bases:
    for i in posFreqMatDict:
        for b in bases:
            try:
                posFreqMatDict[i][b] += 0
            except KeyError:
                posFreqMatDict[i][b] = 0
            
    ## 3. Write out header
    header = ''
    orderedBases = ''
    for e in bases:
        orderedBases += e
        header += e+"\t"
    header = header[:-1]
    print header

    ## 4. write out such that each line is equal to current position
    for i in range(len(posFreqMatDict.keys())):
        line = ''
        for b in orderedBases:
            line += str(posFreqMatDict[i][b])+"\t"
        line = line[:-1]
        print line

def pfmConsensus(posFreqMatDict, verbose=False):
    ''' Returns max of A,C,G,T -- if tie, it randomly selects one'''
    bases = ['A','C','G','T']
    baseDict = {0:'A',1:'C',2:'G',3:'T'}
    cSeq = ''
    numRands = 0
    randIndexes = []
    for i in range(len(posFreqMatDict.keys())):
        options = [posFreqMatDict[i]['A'], posFreqMatDict[i]['C'], posFreqMatDict[i]['G'], posFreqMatDict[i]['T']]
        maxOpt = max(options)
        indexes = [k for k, j in enumerate(options) if j == maxOpt]
        if len(indexes) > 1:
            indexes = [indexes[random.randint(0,len(indexes)-1)]]
            numRands += 1
            randIndexes.append(i)
        cSeq += baseDict[indexes[0]]
        
    sys.stdout.write(cSeq+"\n")
    if verbose:
        sys.stderr.write(str(numRands) + " positions had ties for which a base was randomly selected.\n")
        sys.stderr.write("Positions: " + (", ").join([str(e) for e in randIndexes]) + '\n')
    
    
#def stackMultAln(fastaFile, blocksize=50):
    
