#!/usr/bin/env python2.7

import sys, argparse
from collections import defaultdict
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

DESCRIPTION - combined individual KineticsTools GFFs.



DEFAULTS
MINCOV=25 ## Recommend filtering for >= 25X
MIN10LOGP=20 ## Recommend filtering for >= 20 corresponding to p-valueof 0.01.
                If used a p-value cutoff of 0.01 during modification calling, then 20 should actually already be lowest score in the file.
		However when a identification is attempted, it seems -log10(p) can be < 20.
		Anything labeled only as modified_base is > 20.
		Therefore, all candidate mods are prob >20 before identification, and it becomes less so later.
		I can see that the -log10(p) and identificationQv are NOT the same -- so it is not just simply replaced...
                The p-value prior to identification seems to be from the null hypothesis:
                    that there is no difference in the IPD ratio over this base than the regular model.
                If the p-value here can be interpreted as a regular p-value, then at random one would expect 1 FP per 100 tests.
                    Since it looks at every base in the genome on each strand, there are 2*G tests.
                    Therefore, at random, one could expect 0.01*2*G FPs.
                    In practice, it is hard to tell if this analysis follows that behavior.
                    
MINQV=20  ## Recommend trying >= 20  (1 in 100 expected to be wrong) for higher quality motif analyses (this does not apply to modified_base)



NOTES ON GFF FILE BELOW
Column	Description
Seqid	Reference sequence tag for this observation. Same as refName in the .csv file.
Source	Name of tool -- "kinModCall"
Type	Modification type - a generic tag "modified_base" is used for unidentified bases. For identified bases, m6A, m4C, and m5C are used.
Start	Location of modification.
End	Location of modification.
Score	-10 log (p-value) score for the detection of this event. Analogous to a Phred quality score.
        A value of 20 is the minimum default threshold for this file, and corresponds to a p-value of 0.01.
        A score of 30 corresponds to a p-value of 0.001.
Strand	Native sample strand where kinetics were generated.
        "+" is the strand of the original FASTA and "-" is the reverse complement of the strand. 
	Note that in the .csv file these are marked "0" and "1" respectively.
Phase	Not applicable.
Attributes	Contains extra fields.
                IPDRatio    traditional IPD Ratio
                context     the reference sequence -20bp to +20bp around the modification plus the base at this location as the 21st character,
                            and sequencing coverage of that position.
                            Context is always written in 5' -> 3' orientation of the template strand.
                If the modification type is determined and coverage at that position is at least 10x,
                    there are several additional metrics:
                    frac                an estimate of the fraction of modified reads
                    fracUp/fracLow      the upper and lower 95% confidence intervals 
		    identificationQv    confidence score for the modification type is reported (calculated like the Score).


    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('gff', metavar='gff', nargs='*',
                   type= str, 
                   help='''Path to as many GFFs as need be.''')


parser.add_argument('--modtype', '-m', type=str, default="m5C,m4C,m6A,modified_base", help='''Modtype: can be any single or comma-separated group from m5C, m4C, m6A, modified_base. Default (all): m5C,m4C,m6A,modified_base''')
parser.add_argument('--midbase', '-b', type=str, default="ACGT", help='''This is only used to help further filter modified_base by its midbase.''')
parser.add_argument('--mincov', '-c', type=int, default=25, help='''Minimum coverage to consider modification. Defualt: 25.''')
parser.add_argument('--minlog10p', '-p', type=int, default=20, help='''Minimum -log10(p-value). Default setting is 20. Also, the modification calling process likley alread set the cutoff to 20 (prior to identification). I believe this p-value is w/r/t thenull hypothesis that there is no difference in the IPD ratio over this base than the regular model.''')
parser.add_argument('--minQV', '-q', type=int, default=20, help='''Minimum identification QV to consider. Only pertains to m5C/4mC/6mA. Default=20.''')
parser.add_argument('--minfrac', '-f', type=float, default=float('-inf'), help='''Fraction of reads at this position that are modifided should be at least this (between 0-1). Default=No lower bound. ''')
parser.add_argument('--maxfrac', '-F', type=float, default=float('inf'), help='''Fraction of reads at this position that are modifided should be at most this (between 0-1). Default=No upper bound. ''')
parser.add_argument('--fasta', '-s', action='store_true', default=False, help='''Return FASTA format of the context sequences with annotated headers.''')
parser.add_argument('--table', '-t', action='store_true', default=False, help='''Return tabular format without annotated context sequences.''')
parser.add_argument('--tableseq', '-T', action='store_true', default=False, help='''Return tabular format with annotated context sequences.''')
parser.add_argument('--model', '-M', action='store_true', default=False, help='''Return tabular format model of midbase kmers/counts up to 7mers. This can be used to compare to background frequencies or build PWMs.''')
parser.add_argument('--weighted_model', '-W', action='store_true', default=False, help='''Return tabular format model of midbase kmers/counts up to 7mers. Here counts are a sum of the -10log10p scores associated with them. This can be used to build PWMs. Comparing to random might involve shuffling scores around the file before re-computing. ''')
parser.add_argument('--mididx', '-I', type=int, default=20, help='''This should likely not be changed. It is used to help filter modified_base via its --midbase. ''')

args = parser.parse_args()






class GFF_Entry(object):
    def __init__(self, entry):
        ''' entry       = GFF line unsplit
            attributes  = dict output from GFF._parse_desc()

            The following attributes are anticipated though may not be present:
                coverage
                context
                IPDRatio
                frac
                fracUp
                fracLow
                identificationQv
        '''
        self.entry = entry.split('\t')
        
        self.attributes_dict = self._parse_desc(self.attributes().split(';'))

        self.has_attribute = {}
        
    def _parse_desc(self, desc):
        d = {}
        for e in desc:
            try:
                k,v = e.split('=')
                join = '='
            except: #Although not approp - I have seen some GFF desc terms sep by :
                k,v = e.split(':')
                join = ':'
            d[k] = v
        return d
    
    def seqid(self):
        return self.entry[0]
    def source(self):
        return self.entry[1]
    def type(self):
        return self.entry[2]
    def start(self):
        return int(self.entry[3])
    def end(self):
        return int(self.entry[4])
    def score(self):
        return float(self.entry[5])
    def strand(self):
        return self.entry[6]
    def phase(self):
        return self.entry[7]
    def attributes(self):
        return self.entry[8]
    def get_attribute(self, key):
        if self.in_attributes(key):
            return self.attributes_dict[key]
        return None
    def in_attributes(self, key):
        try:
            return self.has_attribute[key]
        except:
            self.has_attribute[key] = key in self.attributes_dict.keys()
            return self.has_attribute[key]

    def attribute_is(self, key, value):
        return self.get_attribute(key) == value
    
    def context(self):
        '''All expected to have context'''
        return self.get_attribute('context')
    def coverage(self):
        '''All expected to have coverage'''
        return int(self.get_attribute('coverage'))

    def frac(self):
        '''I believe only those with >= 10x coverage expected to have frac/fracUp/fracLow.'''
        try:
            return float(self.get_attribute('frac'))
        except TypeError:
            return None

    def fracUp(self):
        '''I believe only those with >= 10x coverage expected to have frac/fracUp/fracLow.'''
        try:
            return float(self.get_attribute('fracUp'))
        except TypeError:
            return None
        
    def fracLow(self):
        '''I believe only those with >= 10x coverage expected to have frac/fracUp/fracLow.'''
        try:
            return float(self.get_attribute('fracLow'))
        except TypeError:
            return None
        
    def identificationQv(self):
        '''I believe only those with >= 10x coverage and other constraints expected to have this.'''
        try:
            return float(self.get_attribute('identificationQv'))
        except TypeError:
            return None
    
    def length(self):
        return self.end() - self.start() + 1 #self.end() - self.start + 1

    def fasta(self):
        name = '>' + self.seqid() + ':' + str(self.start()) + '-' + str(self.end())
        
        description = [name, 'type:'+str(self.type()), 'strand:'+self.strand(), '-10log10p:'+str(self.score()), 'coverage:'+str(self.coverage()), 'idQv:'+str(self.identificationQv()), 'frac:'+str(self.frac()), 'fracLow:'+str(self.fracLow()), 'fracUp:'+str(self.fracUp())]
        return '\t'.join(description) + '\n' + self.context()

    def table(self, addseq=False):
        description = [self.seqid(), str(self.start()), self.strand(), str(self.type()), str(self.score()), str(self.coverage()), str(self.identificationQv()), str(self.frac()), str(self.fracLow()), str(self.fracUp())]
        if addseq:
            description.append( self.context() )
        return '\t'.join(description)

    def __str__(self):
        return '\t'.join(self.entry).strip()


class MidBaseModel(object):
    def __init__(self, kmers=[1,2,3,4,5,6,7], mid=20):
        self.kmers = kmers
        self.model = {k:{i:defaultdict(int) for i in range(k)} for k in kmers}
        self.mid = mid

    def add(self, gff, score_as_count=False):
        '''gff is a GFF_Entry object'''
        seq = gff.context()
        if score_as_count:
            score = gff.score()
        else:
            score = 1
        for k in self.kmers:
            for i in sorted(self.model[k].keys()):
                left = self.mid - i
                right = self.mid + (k-i)
                self.model[k][i][seq[left:right]] += score

    def __str__(self):
        #return str(self.model)
        out = ''
        for k in self.kmers:
            for i in sorted(self.model[k].keys()):
                for seq in sorted(self.model[k][i].keys()):
                    ## Add 1 to i so it is 1-based
                    out += '\t'.join([str(e) for e in [k, i+1, seq, self.model[k][i][seq]]]) + '\n'
        return out.strip()
    
#### EXECUTE
if len(args.gff) == 0:
    args.gff = ['stdin']
##

def test(gff, args): #mod, mincov, minlog10p, minQv, minfrac, maxfrac
    if gff.type() == "modified_base":
        Qv = args.minQV+1 ## This just prevents the need for a separate return line since this type never has an identificationQv
        midbase = gff.context()[args.mididx].upper() in args.midbase.upper() 
    else:
        Qv = gff.identificationQv()
        midbase = True
    if not gff.in_attributes('frac'):
        frac = 0.0 ## If default of -inf in place, this will be included. It will not be by any other limit.
    else:
        frac = gff.frac()

    return midbase and gff.type() in args.modtype.split(',') and gff.coverage() >= args.mincov and gff.score() >= args.minlog10p and Qv >= args.minQV and frac >= args.minfrac and frac <= args.maxfrac



    
def stdin(fh):
    if fh in ('-', 'stdin') or fh.startswith('<('):
        return True
    return False


if args.model or args.weighted_model:
    midbasemodel = MidBaseModel()


try:
    for fh in args.gff:
        f = sys.stdin if stdin(fh) else open(fh)
        for line in f:
            if not line or line.startswith('#'):
                continue
            gff = GFF_Entry(line)
            if test(gff, args):
                if args.fasta:
                    print gff.fasta()
                elif args.table:
                    print gff.table()
                elif args.tableseq:
                    print gff.table(addseq=True)
                elif args.model:
                    midbasemodel.add(gff)
                elif args.weighted_model:
                    midbasemodel.add(gff,score_as_count=True)
                else:
                    print gff
        f = None if stdin(fh) else f.close()
except IOError:
    pass

if args.model or args.weighted_model:
    print midbasemodel


quit()
########################################################################################################
## INCORPORATE
## --model : just outputs midbase model. Make another script for making models off of bg fastas.
		
##mid=defaultdict(int); 
##allbp=defaultdict(int); 
##alldi=defaultdict(int); 
##alltri=defaultdict(int); 
##all5mer=defaultdict(int); 
##all7mer=defaultdict(int); 
##middimer1=defaultdict(int)
##middimer2=defaultdict(int)
##midtrimer=defaultdict(int)
##mid5mer=defaultdict(int)
##mid7mer=defaultdict(int)
##
##f=sys.stdin; 
##nlines=0
##for line in f:
##	nlines+=1
##	line=line.strip().split();
##	mid[line[0][20]] += 1
##	middimer1[line[0][19:21]] += 1
##	middimer2[line[0][20:22]] += 1
##	midtrimer[line[0][19:22]] += 1
##	mid5mer[line[0][18:23]] += 1
##	mid7mer[line[0][17:24]] += 1
##	for e in line[0]:
##		allbp[e] += 1
##	for i in range(len(line[0])-2+1):
##		alldi[line[0][i:i+2]] += 1
##	for i in range(len(line[0])-3+1):
##		alltri[line[0][i:i+3]] += 1
##	for i in range(len(line[0])-5+1):
##		all5mer[line[0][i:i+5]] += 1
##	for i in range(len(line[0])-7+1):
##		all7mer[line[0][i:i+7]] += 1
##
##
#### if different background model given...
##if len(sys.argv) > 1:
##	def fxn(x):
##		for k in x.keys():
##			x[k] = 1
##		return x
##	
##	#ensure all kmers present have pseudocounts
##	allbp=fxn(allbp); 
##	alldi=fxn(alldi); 
##	alltri=fxn(alltri); 
##	all5mer=fxn(all5mer); 
##	all7mer=fxn(all7mer); 
##	print sys.argv[1]
##	for fa in SeqIO.parse(sys.argv[1], 'fasta'):
##		seq = str(fa.seq)+str(fa.seq.reverse_complement())
##		for e in seq:
##			allbp[e] += 1
##		for i in range(len(seq)-2+1):
##			alldi[seq[i:i+2]] += 1
##		for i in range(len(seq)-3+1):
##			alltri[seq[i:i+3]] += 1
##		for i in range(len(seq)-5+1):
##			all5mer[seq[i:i+5]] += 1
##		for i in range(len(seq)-7+1):
##			all7mer[seq[i:i+7]] += 1
##		
##
##print nlines, "sequences\n"
##
##sumbp = float(sum([e for e in allbp.values()]))
##sumdi = float(sum([e for e in alldi.values()]))
##sumtri = float(sum([e for e in alltri.values()]))
##sum5mer = float(sum([e for e in all5mer.values()]))
##sum7mer = float(sum([e for e in all7mer.values()]))
##summid = float(sum([e for e in mid.values()]))
##summiddi1 = float(sum([e for e in middimer1.values()]))
##summiddi2 = float(sum([e for e in middimer2.values()]))
##summidtri = float(sum([e for e in midtrimer.values()]))
##summid5mer = float(sum([e for e in mid5mer.values()]))
##summid7mer = float(sum([e for e in mid7mer.values()]))
##
##print "Sum of each base in the 41mers"
##for e in allbp.keys():
##	print e, allbp[e], allbp[e]/sumbp
##print
##print "Sum of each base occuring in the middle position"
##d = {}
##for e in mid.keys():
##	d[(mid[e]/summid)/(allbp[e]/sumbp)] = e
##for v in sorted(d.keys(), reverse=True):
##	e=d[v]
##	print e, mid[e], mid[e]/summid,  allbp[e]/sumbp, (mid[e]/summid)/(allbp[e]/sumbp)
##
##print
##print "Sum of dimers w/ midbase as second base"
##d={}
##for e in middimer1.keys():
##	d[(middimer1[e]/summiddi1)/(alldi[e]/sumdi)] = e
##for v in sorted(d.keys(), reverse=True):
##	e=d[v]
##	print e, middimer1[e], middimer1[e]/summiddi1, alldi[e]/sumdi, (middimer1[e]/summiddi1)/(alldi[e]/sumdi)
##
##print
##print "Sum of dimers w/ midbase as first base"
##d={}
##for e in middimer2.keys():
##	d[(middimer2[e]/summiddi2)/(alldi[e]/sumdi)] = e
##for v in sorted(d.keys(), reverse=True):
##	e=d[v]
##	print e, middimer2[e], middimer2[e]/summiddi2, alldi[e]/sumdi, (middimer2[e]/summiddi2)/(alldi[e]/sumdi) 
##
##
##print
##print "Sum of trimers w/ midbase as middle base"
##d={}
##for e in midtrimer.keys():
##	d[(midtrimer[e]/summidtri)/(alltri[e]/sumtri)] = e
##for v in sorted(d.keys(), reverse=True):
##	e=d[v]
##	print e, midtrimer[e], midtrimer[e]/summidtri, alltri[e]/sumtri, (midtrimer[e]/summidtri)/(alltri[e]/sumtri) 
##
##
##print
##print "Sum of 5mers w/ midbase as middle base - enriched > 2-fold"
##d={}
##for e in mid5mer.keys():
##	d[(mid5mer[e]/summid5mer)/(all5mer[e]/sum5mer)] = e
##for v in sorted(d.keys(), reverse=True):
##	e=d[v]
##	if (mid5mer[e]/summid5mer)/(all5mer[e]/sum5mer) >= 2:
##		print e, mid5mer[e], mid5mer[e]/summid5mer, all5mer[e]/sum5mer, (mid5mer[e]/summid5mer)/(all5mer[e]/sum5mer) 
##
##print
##print "Sum of 7mers w/ midbase as middle base - enriched > 5-fold"
##d={}
##for e in mid7mer.keys():
##	d[(mid7mer[e]/summid7mer)/(all7mer[e]/sum7mer)] = e
##for v in sorted(d.keys(), reverse=True):
##	e=d[v]
##	if (mid7mer[e]/summid7mer)/(all7mer[e]/sum7mer) >= 5:
##		print e, mid7mer[e], mid7mer[e]/summid7mer, all7mer[e]/sum7mer, (mid7mer[e]/summid7mer)/(all7mer[e]/sum7mer) 

