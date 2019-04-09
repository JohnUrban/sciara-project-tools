
awk '$1 ~ "^>" {sub(/>/,""); print $1}' ${ASM} > all.names
NODES=~/scratch/taxdb/nodes.dmp; 
NAMES=~/scratch/taxdb/names.dmp; 
taxonomyFromTaxID.py -i all.blastout -nc 1 -tc 2 -no $NODES -na $NAMES > all.tax.out
taxonomy-summarizer.py -i all.tax.out -o all.taxsummary -a all.names 
