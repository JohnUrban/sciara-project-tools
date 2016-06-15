#!/usr/local/bin/python
from collections import defaultdict


def detectEncoding(connection):
    possibilities = set(["Sanger", "Solexa", "Illumina 1.3", "Illumina 1.5", "Illumina 1.8+"])
    notSolexaNor13or15 = "!\"#$%&'()*+,-./0123456789:"
    not13nor15 = ";<=>?"
    not15 = "@AB"
    notSanger = "J"
    notSangerNor18 = "KLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
    f = connection
    f.next(); f.next(); f.next()
    i=1
    while len(possibilities) > 1 and i < 100:
        scoreSet = f.next()
##        print scoreSet
        for score in scoreSet:
##            print score, possibilities
            if score in notSolexaNor13or15:
                try:
                    possibilities.remove("Solexa")
                    possibilities.remove("Illumina 1.3")
                    possibilities.remove("Illumina 1.5")
                except KeyError:
                    pass                   
            elif score in not13nor15:
                try:
                    possibilities.remove("Illumina 1.3")
                    possibilities.remove("Illumina 1.5")
                except KeyError:
                    pass
            elif score in not15:
                try:
                    possibilities.remove("Illumina 1.5")
                except KeyError:
                    pass
            elif score in notSangerNor18:
                try:
                    possibilities.remove("Sanger")
                    possibilities.remove("Illumina 1.8+")
                except KeyError:
                    pass
            elif score in notSanger:
                try:
                    possibilities.remove("Sanger")
                except KeyError:
                    pass
            if len(possibilities) == 1:
                break
##            print possibilities
        f.next(); f.next(); f.next()
        i += 1
    f.seek(0)
    try:
        return list(possibilities)[0]
    except IndexError:
        print "Error.."
                
def rosettaStone(encoding):
    '''Takes in encoding string as put out by detectEncoding()
    returns dictionary that maps symbols to scores.'''
    #define symbols based on encoding scheme
    if encoding == "Illumina 1.8+":
        syms = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
        low, high = 0,41
    elif encoding == "Illumina 1.5":
        syms = "BCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
        low,high = 2,40
    elif encoding == "Illumina 1.3":
        syms = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
        low,high = 0,40
    elif encoding == "Sanger":
        syms = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"
        low, high = 0,40
    elif encoding == "Solexa":
        syms = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
        low, high = -5,40
    symToScore = dict()
    for sym in syms:
        symToScore[sym] = low
        low += 1
    return symToScore

def tableString(symToScore):
    '''Input is output of rosettaStone
    turns it into a string that is printed as a table by sys.stdout.write()'''
    return '\n'.join([str(k) + '\t' + str(symToScore[k]) for k in symToScore])


def readLength(connection):
    read = connection.next()
    read = connection.next()[:-1]
    readLen = len(read)
    connection.seek(0)
    return readLen


def baseQC(connection, symToScore, N=float("inf")):
    ''' profile the base qualities of first N reads
    if N=infinity, it goes through entire file'''
    readLen = readLength(connection)
    scoreHist = {symToScore[sym]:0 for sym in symToScore.keys()}
    posHist = {pos:{symToScore[sym]:0 for sym in symToScore.keys()} for pos in range(1,readLen+1)}   
    f = connection
    i = 1
    while f and i <= N:
        try:
            f.next(); f.next(); f.next()
            symbolString = f.next()[:-1]
            for pos in range(1,readLen+1):
                score = symToScore[symbolString[pos-1]]
                posHist[pos][score] += 1
            i += 1
        except StopIteration: break
    f.seek(0)
    return posHist
    

def baseQCmeans(posHist):
    means = []
    numReads = float(sum(posHist[1].values()))
    for pos in posHist.keys():
        mean = 0
        for score in posHist[pos]:
            mean += posHist[pos][score]*score
        mean = mean/numReads
        means.append(mean)
    return means

def baseQCquantile(posHist, quantile=0.5):
    '''always takes the floor when not an integer -- e.g. 0.25 of 100 and 101 will both give 25th score'''
    quantileScores = []
    numReads = float(sum(posHist[1].values()))
    quantileRank = int(numReads*quantile)
    currRank= 0
    for pos in posHist.keys():
        currRank = 0
        for score in posHist[pos]:
            currRank += posHist[pos][score]
            if currRank > quantileRank:
                break
        quantileScores.append(score-1)
    return quantileScores

def baseQCquantiles(posHist, quantile):
    pass
                

if __name__ == "__main__":
    import sys
    if len(sys.argv) <= 1:
        print "Usage: fqEncoding.py file.fastq"
        quit()
    connection = open(sys.argv[1], "r")
    print detectEncoding(connection)
