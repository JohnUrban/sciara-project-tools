## run after blast analysis complete - 20181216

cat blast/* > canu_windows.blastout  ## This file can be used for a window-specific blob analysis (or just the first 3 cols) 
## -- need to get cov for each window, be sure to normalize by window length so the edge cases that are not 1000 bp do not cause too much noise -- Or only look at windows that are 1000 bp (i.e. skip the last window of each contig)

awk 'OFS="\t" {gsub(/\ /,"_"); gsub(/__/,""); sub(/:/,"\t"); sub(/-/,"\t"); print $1,$4,$5}' canu_windows.blastout | sort -k1,1 -k2,2n -k3,3nr > canu_blastn_NT_qseqTaxBit.txt

echo "#chr=qseqid start=qstart end=qend name=gi;taxID;bitscore;pident;length;mismatch;gapopen;sstart;send;evalue;sscinames;sskingdoms;stitle" > canu_blastn_NT_annotation.bed
awk 'OFS="\t" {gsub(/\ /,"_"); gsub(/__/,""); sub(/:/,"\t"); sub(/-/,"\t"); print $1,$2+$11-1,$2+$12,$6";"$4";"$5";"$7";"$8";"$9";"$10";"$13"-"$14";"$15";"$16";"$17";"$18}' canu_windows.blastout | sortBed -i - >> canu_blastn_NT_annotation.bed


## some things for blob excluding bradysia hits
grep -v -E -i 'Bradysia|coprophila|DASH' canu_windows.blastout | awk 'OFS="\t" {gsub(/\ /,"_"); gsub(/__/,""); sub(/:/,"\t"); sub(/-/,"\t"); print $1,$4,$5}' | sort -k1,1 -k2,2n -k3,3nr > canu_blastn_NT_qseqTaxBit_excluding-Bradysia-coprophila-DASH.txt
grep -v DASH canu_windows.blastout | awk 'OFS="\t" {gsub(/\ /,"_"); gsub(/__/,""); sub(/:/,"\t"); sub(/-/,"\t"); print $1,$4,$5}' | sort -k1,1 -k2,2n -k3,3nr > canu_blastn_NT_qseqTaxBit_excluding-DASH.txt

####
exit
#### Looking at DASH prevalence here...
grep DASH canu_blastn_NT_annotation.bed | mergeBed -i - | awk '{s+=$3-$2}END{print s,NR,s/NR}'
10517746 27484 382.686

grep DASH canu_blastn_NT_annotation.bed | mergeBed -i - | intersectBed -u -a - -b <(grep -v DASH canu_blastn_NT_annotation.bed) | awk '{s+=$3-$2}END{print s,NR,s/NR}'
2799499 2509 1115.78

grep DASH canu_blastn_NT_annotation.bed | mergeBed -i - | intersectBed -v -a - -b <(grep -v DASH canu_blastn_NT_annotation.bed) | awk '{s+=$3-$2}END{print s,NR,s/NR}'
7718247 24975 309.039

grep DASH canu_blastn_NT_annotation.bed | mergeBed -i - | intersectBed -u -a <(grep -v DASH canu_blastn_NT_annotation.bed) -b - | awk '{s+=$3-$2}END{print s,NR,s/NR}'
3245798 13643 237.909

