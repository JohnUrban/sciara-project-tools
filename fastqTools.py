import random
import sys
import os
import numpy as np
from collections import defaultdict, deque


def bernoulli(probability):
    return np.random.binomial(1, float(probability))



def fastqToRaw(fqFile):
    ''' fastq structure is 4 lines per unit
        lines 1-4 == @name, sequence, +name, quality scores
        This extracts just the sequences (line 2s)
        it prints each so as to not store in memory
        Thus, this needs to be used from command line and redirected'''
    ctr = 0
    for line in open(fqFile, 'r'):
        ctr += 1
        if ctr%4 == 2:
            print line[:-1]


def numReads(fqFile):
    with open(fqFile, 'r') as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return int((i+1)/4)


def extractRead(fqFile, readNum):
    '''the readNumTh is the read that occurs at line readNum*4'''
    ctr = 0
    if readNum == 0:
        totalNumRead = numReads(fqFile)
        readNum = random.randint(1,totalNumRead)
    f = open(fqFile, 'r')
    try:
        while f:
            ctr += 1
            if ctr == readNum:
                read = f.next() + f.next() + f.next() + f.next()
                sys.stdout.write(read)
                break
            else:
                for i in range(4): f.next()
    except StopIteration:
        pass


def shuffleFqFile(fqFile, outputFile=None, givenTotalNumReadsInFile=None):
    ''' 1 M reads will take up >= 3*8 MB = 24 MB memory
        10 M will take uup 240 MB
        100 M will take up 2.4 GB
        200 M will take up >= 4.8 GB
        Thus, this is probably only useful on files with < 100 M reads.
        A wrapper can be made that divides larger files up, shuffles smaller files, and randomly concats them back.'''
    if not givenTotalNumReadsInFile:
        totalNumReads = numReads(fqFile)
    else:
        totalNumReads = givenTotalNumReadsInFile
    readIndexes = deque(range(1, totalNumReads+1))
        ## need this to be deque instead of list b/c lists are fast only for appending/popping items off end
        ## deques are fast at popping them off both sides: beginning and end (could also just choose to pop off end with list)
    random.shuffle(readIndexes)
    ## index file: open index file and read in forming dictionary with k:v == readNum:startInMem
    intIndex = internalIndex(fqFile)
    ## go through shuffled list (pop stuff off to shrink it), go to read location, extract read, and print read
    f = open(fqFile, 'r')
    if not outputFile:
        while readIndexes:
            f.seek(intIndex[readIndexes[0]])
            fastqEntry = f.next() + f.next() + f.next() + f.next()
            sys.stdout.write(fastqEntry)
            readIndexes.popleft()
    else: ## outputfile specified
        out = open(outputFile, 'w')
        while readIndexes:        
            f.seek(intIndex[readIndexes[0]])
            fastqEntry = f.next() + f.next() + f.next() + f.next()
            out.write(fastqEntry)
            readIndexes.popleft()      



def internalIndex(fqFile):
    '''This is the same as the "index" function only it is made for functions to use internally.
        Whereas index has O(constant) space req, this has O(n) space req.
        This is b/c "index" just goes through a file and writes to another file as it goes.
        And "internalIndex" stores a dictionary with n k,v pairs where n is number of reads.
        The k,v pairs are readIndex and byte it starts at (where a read's index is k if it is the kth read in file)
        "internalIndex" returns the dict
        Another difference is that index returns positions of actual sequence whereas this returns positions of read name
        '''
    f = open(fqFile, 'rb')
    f.seek(0)
    ## initial
    ctr = 0
    fqEntry = f.next() + f.next() + f.next() + f.next()
    size = sys.getsizeof(fqEntry) - 37
    sumSize = 1 
    intIndex = defaultdict(int)
    while fqEntry:
        ctr += 1
        intIndex[ctr] = sumSize ## the sum size when up to a 'line 1' is exact stating pos of a 'line2' b/c start is 0 index and bytes are 1-index
        try:
            fqEntry = f.next() + f.next() + f.next() + f.next()
        except StopIteration:
            break
        size = sys.getsizeof(fqEntry) - 37
        sumSize += size
    return intIndex
        

def downSampleReadsWithoutReplacement(fqFile, probability, outputFile=None):
    f = open(fqFile, 'r')
    if not outputFile:
        try:
            while f:
                read = f.next() + f.next() + f.next() + f.next()
                if bernoulli(probability):
                    sys.stdout.write(read)
                    ## this is used instead of "print" b/c it doesn't add extra newlines
                    ## print read[:-1] seemed to work too, but I feel more comfortable with this
        except StopIteration:
            pass
    else:
        out = open(outputFile, 'w')
        try:
            while f:
                read = f.next() + f.next() + f.next() + f.next()
                if bernoulli(probability):
                    out.write(read)
                    ## this is used instead of "print" b/c it doesn't add extra newlines
                    ## print read[:-1] seemed to work too, but I feel more comfortable with this
        except StopIteration:
            pass        

def downSampleReadsWithReplacement(fqFile, numReadsToSample, outputFile=None):
    numReadsInFile = numReads(fqFile)
    getThisRead = defaultdict(int)
    for i in range(numReadsToSample):
        getThisRead[random.randint(1,numReadsInFile)] += 1
    readNum = 0
    f = open(fqFile, 'r')
    if not outputFile: ## stdout
        try:
            while f:
                readNum += 1
                read = f.next() + f.next() + f.next() + f.next()
                while getThisRead[readNum] > 0:
                    sys.stdout.write(read)
                    getThisRead[readNum] -= 1
        except StopIteration:
            pass
    else: ## write to file
        out = open(outputFile, 'w')
        try:
            while f:
                readNum += 1
                read = f.next() + f.next() + f.next() + f.next()
                while getThisRead[readNum] > 0:
                    out.write(read)
                    getThisRead[readNum] -= 1
        except StopIteration:
            pass                              
        
    

def estDistBtwPatInGenomeFromString(string, pattern):
    ## estimateDistanceBetweenPatternsInGenomeFromFq
    ''' This will go through a strinf and estimate the mean and median distance
        between instances of a given pattern (e.g. restriction site) in the genome
        This fxn goes line by line (i.e. does not merge reads)
        output: numDistances, meanDist, medianDist'''
    patLen = len(pattern)
    dist = 0
    distances = []
    for i in range(len(string) - patLen):
        dist += 1
        if string[i:i+patLen] == pattern:
            distances += [dist]
            dist = 0
    ## analyze distances
    return analyzeMeansAndMedians(distances)

def estDistBtwPatInGenomeFromFq(fqFile, pattern):
    ## estimateDistanceBetweenPatternsInGenomeFromFq
    ''' This will go through fastq file and estimate the mean and median distance
        between instances of a given pattern (e.g. restriction site) in the genome
        This fxn goes line by line (i.e. does not merge reads)
        output: numDistances, meanDist, medianDist'''
    ctr = 0
    patLen = len(pattern)
    dist = 0
    distances = []
    for read in open(fqFile, 'r'):
        ctr += 1
        if ctr%4 == 2:
            read = read[:-1] ## cutting of \n at end of line ("read")
            for i in range(len(read) - patLen):
                dist += 1
                if read[i:i+patLen] == pattern:
                    distances += [dist]
                    dist = 0
    ## analyze distances
    return analyzeMeansAndMedians(distances)
            

def analyzeMeansAndMedians(distances):
    ## analyze distances
    distances = sorted(distances)
    sumDistances = float(sum(distances))
    numDistances = len(distances)
    # calc mean
    try:
        mean = sumDistances/numDistances
    except ZeroDivisionError:
        mean = None
    ## calc median
    try:
        if numDistances%2 == 0: ## even
            median = (distances[numDistances//2 - 1] + distances[(numDistances//2)])/2.0
        else: ## odd
            median = distances[(numDistances//2)]
    except IndexError:
        median = None
    if numDistances < 2:
        mean = None
    return numDistances, mean, median

def estDistBtwPatBySamplingInfinitelyFromFqFile(fqFile, pattern, numReadsToSample, replacement=True, givenTotalNumReadsInFile=None, intIndex=None):
    ''' This will sample with or without replacement
        It takes the approach of just picking a random number 1 at a time
        and sampling the corresponding read while keeping track of dist in the process'''
    ## index file: open index file and read in forming dictionary with k:v == readNum:startInMem
    if not intIndex:
        print "indexing..."
        intIndex = internalIndex(fqFile)
    if not givenTotalNumReadsInFile:
        print "getting numReads..."
        totalNumReads = len(intIndex.keys()) ##numReads(fqFile)
    else:
        totalNumReads = givenTotalNumReadsInFile
    patLen = len(pattern)

    ## if sample w/o replacement, N can be > numReads BUT you have to sample all reads before sampling 1 again
    f = open(fqFile, 'r')
    print "executing..."
    if not replacement:
##        print numReadsToSample, totalNumReads
        numLoops = numReadsToSample/float(totalNumReads)
        numReadsSampled = 0
        distances = []
        while numLoops > 0:
##            print numLoops
            numLoops -= 1
            readIndexes = deque(range(1, totalNumReads+1))
            random.shuffle(readIndexes)
            dist = 0
            while readIndexes and numReadsSampled < numReadsToSample:
                numReadsSampled += 1
                f.seek(intIndex[readIndexes[0]]) ## go to this read
                readIndexes.popleft() ## remove this index
                f.next() ## read is on next line
                read = f.next()[:-1] ## take entire read but not line break
                for i in range(len(read) - patLen):
                    dist += 1
                    if read[i:i+patLen] == pattern:
                        distances += [dist]
                        dist = 0
    ## w/ rep
    else: # WITH REPLACEMENT
        dist = 0
        distances = []
        numReadsSampled = 0
        while numReadsToSample > 0:
            numReadsToSample -= 1
##            numReadsSampled += 1
            readIndex = random.randint(1,totalNumReads)
            f.seek(intIndex[readIndex])
            f.next()
            read = f.next()[:-1]
            for i in range(len(read) - patLen):
                dist += 1
                if read[i:i+patLen] == pattern:
                    distances += [dist]
                    dist = 0

    ## analyze distances
    return analyzeMeansAndMedians(distances)




    
def estNumSitesPer(fqFile, pattern, perDist=100e+3):
    ## related to estimateDistanceBetweenPatternsInGenomeFromFq
    ''' '''
    ctr = 0
    patLen = len(pattern)
    numPatPerDist = 0
    counts = []
    dist = 0
    for read in open(fqFile, 'r'):
        ctr += 1
        if ctr%4 == 2:
            read = read[:-1] ## cutting of \n at end of line ("read")
##            print read
            for i in range(len(read) - patLen):
                dist += 1
##                print dist, perDist
                if read[i:i+patLen] == pattern:
                    numPatPerDist += 1
                if dist == perDist:
                    counts += [numPatPerDist]
                    numPatPerDist = 0
                    dist = 0
    ## analyze numCounts
    return analyzeMeansAndMedians(counts)


def estDensityBySamplingInfinitelyFromFqFile(fqFile, pattern, numReadsToSample, perDist=100e+3, replacement=True, givenTotalNumReadsInFile=None):
    ''' This will sample with or without replacement
        It takes the approach of just picking a random number 1 at a time
        and sampling the corresponding read while keeping track of dist in the process'''
    if not givenTotalNumReadsInFile:
        totalNumReads = numReads(fqFile)
    else:
        totalNumReads = givenTotalNumReadsInFile
    patLen = len(pattern)
    ## index file: open index file and read in forming dictionary with k:v == readNum:startInMem
    intIndex = internalIndex(fqFile)
    
    ## if sample w/o replacement, N can be > numReads BUT you have to sample all reads before sampling 1 again
    f = open(fqFile, 'r')
    numReadsSampled = 0
    counts = []
    numPatPerDist = 0
    dist = 0
##    print "HI"
    if not replacement:
##        print numReadsToSample, totalNumReads
        numLoops = numReadsToSample/float(totalNumReads)
        while numLoops > 0:
##            print numLoops
            numLoops -= 1
            readIndexes = deque(range(1, totalNumReads+1))
            random.shuffle(readIndexes)
            dist = 0
            while readIndexes and numReadsSampled < numReadsToSample:
                numReadsSampled += 1
                f.seek(intIndex[readIndexes[0]]) ## go to this read
                readIndexes.popleft() ## remove this index
                f.next() ## read is on next line
                read = f.next()[:-1] ## take entire read but not line break
                for i in range(len(read) - patLen):
                    dist += 1
                    if read[i:i+patLen] == pattern:
                        numPatPerDist += 1
                    if dist == perDist:
                        counts += [numPatPerDist]
                        numPatPerDist = 0
                        dist = 0
    else: # WITH REPLACEMENT
        while numReadsToSample > 0:
            numReadsToSample -= 1
##            numReadsSampled += 1
            readIndex = random.randint(1,totalNumReads)
            f.seek(intIndex[readIndex])
            f.next()
            read = f.next()[:-1]
            for i in range(len(read) - patLen):
                dist += 1
                if read[i:i+patLen] == pattern:
                    numPatPerDist += 1
                if dist == perDist:
                    counts += [numPatPerDist]
                    numPatPerDist = 0
                    dist = 0

    ## analyze counts
    return analyzeMeansAndMedians(counts)



def GCcontent(): pass




def generateString(length, alphabet=['A','C','G','T']):
    string = ''
    alphLen = len(alphabet)
    for i in range(length):
        string += alphabet[random.randint(0,alphLen-1)]
    return string

def complement(DNAstring):
    DNAstring = DNAstring.upper()
    compString = ''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in DNAstring:
        compString = compString + complement[base]
    return compString

def reverseComplement(DNAstring):
    return complement(DNAstring[-1::-1])



def index(fqFile, indexFile):
    ## emptry string in python costs 37 bytes
    ## each character thereafter costs 1 byte
    ## This makes a file of the starting byte of each seq line
    ## for now its stored as a text file -- reucing the number of lines to iterate over by 4fold
    ## either I will need to do reading in tricks to balance mem and speed (.e. reading in multiple starts at once)
    ## or store the index as a pickle or json file to load back when using
    f = open(fqFile, 'rb')
    ## get readLength
    ctr = 0
    for line in f:
        ctr += 1
        if ctr%4 == 2:
            readLen = len(line[:-1])
            break
    f.seek(0)
    index = open(indexFile,'w')
    index.write("readLen="+str(readLen)+"\n")

    ## initial
    ctr = 0
    line = f.next()
    size = sys.getsizeof(line) - 37
    sumSize = size
    while line:
        ctr += 1
        if ctr%4 == 1: ## then its the next byte
            index.write(str(sumSize) + "\n") ## the sum size when up to a 'line 1' is exact stating pos of a 'line2' b/c start is 0 index and bytes are 1-index
        try:
            line = f.next()
        except StopIteration:
            break
        size = sys.getsizeof(line) - 37
        sumSize += size
##        print size, sumSize
        
        
def readIndex(fqFile, indexFile):
    f = open(fqFile, 'r')
    g = open(indexFile, 'r')
    readLen = int(g.next().split("=")[1])
##    print readLen
    for line in g:
        f.seek(int(line[:-1]))
##        read = f.read(readLen)
##        print "line", line, "read", read
        yield f.read(readLen)
               



def estDistBtwPatInGenomeFromIndex1(fqFile, indexFile, pattern):
    ## estimateDistanceBetweenPatternsInGenomeFromFq
    ''' This will go through fastq file and estimate the mean and median distance
        between instances of a given pattern (e.g. restriction site) in the genome
        This fxn goes line by line (i.e. does not merge reads)
        output: numDistances, meanDist, medianDist'''
    ctr = 0
    patLen = len(pattern)
    dist = 0
    distances = []
    for read in readIndex(fqFile, indexFile):
##        ctr += 1
##        print read
        for i in range(len(read) - patLen):
            dist += 1
            if read[i:i+patLen] == pattern:
                distances += [dist]
                dist = 0
    ## analyze distances
    distances = sorted(distances)
    sumDistances = float(sum(distances))
    numDistances = len(distances)
    # calc mean
    try:
        mean = sumDistances/numDistances
    except ZeroDivisionError:
        mean = None
    ## calc median
    try:
        if numDistances%2 == 0: ## even
            median = (distances[numDistances//2 - 1] + distances[(numDistances//2)])/2.0
        else: ## odd
            median = distances[(numDistances//2)]
    except IndexError:
        median = None
    if numDistances < 2:
        mean = None
    print numDistances, mean, median
    return numDistances, mean, median



def estNumSitesPerFromIndex(fqFile, indexFile, pattern, perDist=100e+3):
    ## related to estimateDistanceBetweenPatternsInGenomeFromFq
    ''' This will go through fastq file and estimate the mean and median distance
        between instances of a given pattern (e.g. restriction site) in the genome
        This fxn goes line by line (i.e. does not merge reads)
        output: numDistances, meanDist, medianDist'''
##    ctr = 0
    patLen = len(pattern)
    numPatPerDist = 0
    counts = []
    dist = 0
    for read in readIndex(fqFile, indexFile):
##        ctr += 1
##        print read
        for i in range(len(read) - patLen):
            dist += 1
##                print dist, perDist
            if read[i:i+patLen] == pattern:
                numPatPerDist += 1
            if dist == perDist:
                counts += [numPatPerDist]
                numPatPerDist = 0
                dist = 0
    ## analyze distances
    counts = sorted(counts)
    sumCounts = float(sum(counts))
    numCounts = len(counts)
    # calc mean
    try:
        mean = sumCounts/numCounts
    except ZeroDivisionError:
        mean = None
    ## calc median
    try:
        if numCounts%2 == 0: ## even
            median = (counts[numCounts//2 - 1] + counts[(numCounts//2)])/2.0
        else: ## odd
            median = counts[(numCounts//2)]
    except IndexError:
        median = None
    if numCounts < 2:
        mean = None
    print numCounts, mean, median
    return numCounts, mean, median



### experimental
def readIndex2(indexFile):
    g = open(indexFile, 'r')
    readLen = int(g.next().split("=")[1])
    readPos = []
    for line in g:
        readPos.append(int(line[:-1]))
    return readLen, readPos
                    
def estDistBtwPatInGenomeFromIndex2(fqFile, indexFile, pattern):
    ## estimateDistanceBetweenPatternsInGenomeFromFq
    ''' This will go through fastq file and estimate the mean and median distance
        between instances of a given pattern (e.g. restriction site) in the genome
        This fxn goes line by line (i.e. does not merge reads)
        output: numDistances, meanDist, medianDist'''
    ctr = 0
    patLen = len(pattern)
    dist = 0
    distances = []
    readLen, readPos = readIndex2(indexFile)
##    print sys.getsizeof(readPos)
    f = open(fqFile, 'r')
    for pos in readPos:
        f.seek(pos)
        read = f.read(readLen)
##        ctr += 1
##        print read
        for i in range(len(read) - patLen):
            dist += 1
            if read[i:i+patLen] == pattern:
                distances += [dist]
                dist = 0
    ## analyze distances
    distances = sorted(distances)
    sumDistances = float(sum(distances))
    numDistances = len(distances)
    # calc mean
    try:
        mean = sumDistances/numDistances
    except ZeroDivisionError:
        mean = None
    ## calc median
    try:
        if numDistances%2 == 0: ## even
            median = (distances[numDistances//2 - 1] + distances[(numDistances//2)])/2.0
        else: ## odd
            median = distances[(numDistances//2)]
    except IndexError:
        median = None
    if numDistances < 2:
        mean = None
    print numDistances, mean, median
    return numDistances, mean, median



### need to run simulation that takes both fwd and reverse reads (in ~equiv amounts)
        ## and see whether searching for a non-palindrome density is better estimated by
        ## 1. searching only the sequence in the reads
        ## or
        ## 2. both the seq and its rev comp
        ## since the reads are both fwd and rev -- it might be just for the seq


def generatePracticeFastqFile(genome, numReads, readLen, outputFqFileName, includeRevComp=False):
##    genome = generateString(genomeLen)
    genomeLen = len(genome)
    if includeRevComp:
        revGenome = reverseComplement(genome) ## assuming the genome length in sim will be under 1 Gbp (which would require ~1GB+37byte mem and ~2GB for btoh)
    reads = []
    f = open(outputFqFileName, 'w')
    for i in xrange(numReads):
        fwd = True
        if includeRevComp:
            fwd = random.randint(0,1)
        index = random.randint(0, genomeLen - 1 - readLen)
        if fwd:
##            print "fwd"
            read = genome[index:index+readLen]
        else: ## rev
##            print "rev"
            read = revGenome[index:index+readLen]
        f.writelines("@read_#_" + str(i+1) + "\n")
        f.writelines(read + "\n")
        f.writelines("+" + "\n")
        f.writelines("I"*readLen + "\n") ## illumina 1.8 encoding where I = qual scor eof 41 (best)
            

def simulateDistanceFromGenomicReads(genome, numReads, readLen, pattern, includeRevComp=False):
    ''' This samples both fwd and rev comp reads to simulate real data.
    It assumes there are ~50% of both types of read.
    It will either report mean and median distances for just between fwd patterns
    or it will look for both fwd and revComp patterns to take distances (includeRevComp=True)'''
    # Begin
    outputFqFileName = "distance-simulation-" + pattern + "-" + str(random.randint(1000,9999)) + ".fq"
    generatePracticeFastqFile(genome, numReads, readLen, outputFqFileName, includeRevComp)
    numDist, mean, median = estDistBtwPatInGenomeFromFq(outputFqFileName, pattern)
    os.system("rm " + outputFqFileName)
    return mean, median

def multiSimDistFromGenomicReads(numSim, genomeLen, numReads, readLen, pattern, includeRevComp=False):
    means = []
    medians = []
    genome = generateString(genomeLen)
    numDist, actualMean, actualMedian = estDistBtwPatInGenomeFromString(genome, pattern)
    for i in range(numSim):
        mean, median = simulateDistanceFromGenomicReads(genome, numReads, readLen, pattern, includeRevComp)
        means.append(mean)
        medians.append(median)
    estMean = 1.0*sum(means)/len(means)
    estMedian = 1.0*sum(medians)/len(medians)
##    medians.sort()
##    numMedians = len(medians)
##    if numMedians > 1:
##        if numMedians%2 == 0:
##            estMedian = (medians[numMedians//2 - 1] + medians[numMedians//2])/2.0
##        else: #odd
##            estMedian = medians[numMedians//2]
##    elif numMedians == 1:
##        estMedian = medians[0]
    return actualMean, actualMedian, estMean, estMedian
        

def multiSimDistFromSingleReadsFile(numSim, genomeLen, numReadsFromGenome, readLen, pattern, numReadsToSample, includeRevComp=False, withReplacement=False):
    means = []
    medians = []
    genome = generateString(genomeLen)
    numDist, actualMean, actualMedian = estDistBtwPatInGenomeFromString(genome, pattern)
    ## single reads file
    outputFqFileName = "distance-simulation-" + pattern + "-" + str(random.randint(1000,9999)) + ".fq"
    generatePracticeFastqFile(genome, numReadsFromGenome, readLen, outputFqFileName, includeRevComp)
    sampledReads = "sampledReads-" + pattern + "-" + str(random.randint(1000,9999)) + ".fq"
    shuffledReads = "shuffledSampledReads-" + pattern + "-" + str(random.randint(1000,9999)) + ".fq"
    if withReplacement:
        for i in range(numSim):
            downSampleReadsWithReplacement(outputFqFileName, numReadsToSample, outputFile=sampledReads)
            shuffleFqFile(sampledReads, outputFile=shuffledReads, givenTotalNumReadsInFile=numReadsToSample)
            numDist, mean, median = estDistBtwPatInGenomeFromFq(shuffledReads, pattern)
            means.append(mean)
            medians.append(median)
        os.system("rm " + shuffledReads)    
    else: ##w/o replacement
        if numReadsToSample > numReadsFromGenome:
            return "When sampling w/o replacement, numReadsToSample need be <= numReadsFromGenome"
        probability = numReadsToSample/float(numReadsFromGenome)
##        print probability
        for i in range(numSim):
            downSampleReadsWithoutReplacement(outputFqFileName, probability, outputFile=sampledReads)
            numDist, mean, median = estDistBtwPatInGenomeFromFq(sampledReads, pattern)
            means.append(mean)
            medians.append(median)
    estMean = 1.0*sum(means)/len(means)
    estMedian = 1.0*sum(medians)/len(medians)
    os.system("rm " + outputFqFileName)
    os.system("rm " + sampledReads)
    return actualMean, actualMedian, estMean, estMedian
                                                                    

### The sampling WITH replacement needs to shuffle reads after it samples
### Otherwise, multiple identical reads are stacked in a row and it messes up the estimates


### I should make a sampler that ONLY shuffles reads from single file (rather than bootstrapping)



"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"

##def detectEncoding(fqFile):
##    possibilities = ["Sanger", "Solexa", "Illumina 1.3", "Illumina 1.5", "Illumina 1.8+"]
##    notSolexaNor1.3or1.5 = "!\"#$%&'()*+,-./0123456789:"
##    not1.3nor1.5 = ";<=>?"
##    not1.5 = "@AB"
##    notSanger = "J"
##    notSangerNor1.8 = "KLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
##    f = open(fqFile, "r")
##    f.next(); f.next(); f.next()
##    print f.next()
##    quit()
####    while len(possibilities) > 1:
        
