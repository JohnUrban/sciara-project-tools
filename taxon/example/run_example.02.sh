
## Replace these with /path/to/each if that is not where they are on your system
NODES=~/data/ncbi/taxonomy/20160623/nodes.dmp
NAMES=~/data/ncbi/taxonomy/20160623/names.dmp

if [ ! -f $NODES ]; then echo "NODES file does not exist"; exit; fi
if [ ! -f $NAMES ]; then echo "NAMES file does not exist"; exit; fi

../taxonomyFromTaxID.py -i example.02.blast -nc 1 -tc 2 -no $NODES -na $NAMES > example.02.taxonomy.out
../taxonomy-summarizer.py -i example.02.taxonomy.out -o example.02.taxsum

