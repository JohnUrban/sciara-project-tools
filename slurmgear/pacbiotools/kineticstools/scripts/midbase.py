import sys; 
from collections import defaultdict as dd; 
from Bio import SeqIO
		
mid=dd(int); 
allbp=dd(int); 
alldi=dd(int); 
alltri=dd(int); 
all5mer=dd(int); 
all7mer=dd(int); 
middimer1=dd(int)
middimer2=dd(int)
midtrimer=dd(int)
mid5mer=dd(int)
mid7mer=dd(int)

f=sys.stdin; 
nlines=0
for line in f:
	nlines+=1
	line=line.strip().split();
	mid[line[0][20]] += 1
	middimer1[line[0][19:21]] += 1
	middimer2[line[0][20:22]] += 1
	midtrimer[line[0][19:22]] += 1
	mid5mer[line[0][18:23]] += 1
	mid7mer[line[0][17:24]] += 1
	for e in line[0]:
		allbp[e] += 1
	for i in range(len(line[0])-2+1):
		alldi[line[0][i:i+2]] += 1
	for i in range(len(line[0])-3+1):
		alltri[line[0][i:i+3]] += 1
	for i in range(len(line[0])-5+1):
		all5mer[line[0][i:i+5]] += 1
	for i in range(len(line[0])-7+1):
		all7mer[line[0][i:i+7]] += 1


if len(sys.argv) > 1:
	allbp=dd(int); 
	alldi=dd(int); 
	alltri=dd(int); 
	all5mer=dd(int); 
	all7mer=dd(int); 
	for fa in SeqIO.parse(sys.argv[1], 'fasta'):
		seq = str(fa.seq)+str(fa.seq.reverse_complement())
		for e in seq:
			allbp[e] += 1
		for i in range(len(seq)-2+1):
			alldi[seq[i:i+2]] += 1
		for i in range(len(seq)-3+1):
			alltri[seq[i:i+3]] += 1
		for i in range(len(seq)-5+1):
			all5mer[seq[i:i+5]] += 1
		for i in range(len(seq)-7+1):
			all7mer[seq[i:i+7]] += 1
		

print nlines, "sequences\n"

sumbp = float(sum([e for e in allbp.values()]))
sumdi = float(sum([e for e in alldi.values()]))
sumtri = float(sum([e for e in alltri.values()]))
sum5mer = float(sum([e for e in all5mer.values()]))
sum7mer = float(sum([e for e in all7mer.values()]))
summid = float(sum([e for e in mid.values()]))
summiddi1 = float(sum([e for e in middimer1.values()]))
summiddi2 = float(sum([e for e in middimer2.values()]))
summidtri = float(sum([e for e in midtrimer.values()]))
summid5mer = float(sum([e for e in mid5mer.values()]))
summid7mer = float(sum([e for e in mid7mer.values()]))

print "Sum of each base in the 41mers"
for e in allbp.keys():
	print e, allbp[e], allbp[e]/sumbp
print
print "Sum of each base occuring in the middle position"
d = {}
for e in mid.keys():
	d[(mid[e]/summid)/(allbp[e]/sumbp)] = e
for v in sorted(d.keys(), reverse=True):
	e=d[v]
	print e, mid[e], mid[e]/summid,  allbp[e]/sumbp, (mid[e]/summid)/(allbp[e]/sumbp)

print
print "Sum of dimers w/ midbase as second base"
d={}
for e in middimer1.keys():
	d[(middimer1[e]/summiddi1)/(alldi[e]/sumdi)] = e
for v in sorted(d.keys(), reverse=True):
	e=d[v]
	print e, middimer1[e], middimer1[e]/summiddi1, alldi[e]/sumdi, (middimer1[e]/summiddi1)/(alldi[e]/sumdi)

print
print "Sum of dimers w/ midbase as first base"
d={}
for e in middimer2.keys():
	d[(middimer2[e]/summiddi2)/(alldi[e]/sumdi)] = e
for v in sorted(d.keys(), reverse=True):
	e=d[v]
	print e, middimer2[e], middimer2[e]/summiddi2, alldi[e]/sumdi, (middimer2[e]/summiddi2)/(alldi[e]/sumdi) 


print
print "Sum of trimers w/ midbase as middle base"
d={}
for e in midtrimer.keys():
	d[(midtrimer[e]/summidtri)/(alltri[e]/sumtri)] = e
for v in sorted(d.keys(), reverse=True):
	e=d[v]
	print e, midtrimer[e], midtrimer[e]/summidtri, alltri[e]/sumtri, (midtrimer[e]/summidtri)/(alltri[e]/sumtri) 


print
print "Sum of 5mers w/ midbase as middle base - enriched > 2-fold"
d={}
for e in mid5mer.keys():
	d[(mid5mer[e]/summid5mer)/(all5mer[e]/sum5mer)] = e
for v in sorted(d.keys(), reverse=True):
	e=d[v]
	if (mid5mer[e]/summid5mer)/(all5mer[e]/sum5mer) >= 2:
		print e, mid5mer[e], mid5mer[e]/summid5mer, all5mer[e]/sum5mer, (mid5mer[e]/summid5mer)/(all5mer[e]/sum5mer) 

print
print "Sum of 7mers w/ midbase as middle base - enriched > 5-fold"
d={}
for e in mid7mer.keys():
	d[(mid7mer[e]/summid7mer)/(all7mer[e]/sum7mer)] = e
for v in sorted(d.keys(), reverse=True):
	e=d[v]
	if (mid7mer[e]/summid7mer)/(all7mer[e]/sum7mer) >= 5:
		print e, mid7mer[e], mid7mer[e]/summid7mer, all7mer[e]/sum7mer, (mid7mer[e]/summid7mer)/(all7mer[e]/sum7mer) 

