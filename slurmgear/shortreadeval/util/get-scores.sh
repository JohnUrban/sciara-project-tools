#!/bin/bash

export PATH=/gpfs/scratch/jurban/male-ilmn/long_read_evals/scripts/:$PATH

## after all renaming and such that was necessary to have thngs sort in same orders etc

source score-fxns.sh

##Score top dir
echo Top
scoredir

##Score other dirs
for dir in canu_quiver_as_is original_recipe; do
 cd $dir
 echo $dir
 scoredir
 cd ../
done

unitescores quiver0/combined-scores.quiver0.txt quiver1/combined-scores.quiver1.txt quiver2/combined-scores.quiver2.txt quiver3/combined-scores.quiver3.txt quiver4/combined-scores.quiver4.txt quiver5/combined-scores.quiver5.txt quiver6/combined-scores.quiver6.txt quiver7/combined-scores.quiver7.txt

cd original_recipe
unitescores quiver1/combined-scores.quiver1.txt quiver2/combined-scores.quiver2.txt quiver3/combined-scores.quiver3.txt quiver4/combined-scores.quiver4.txt quiver5/combined-scores.quiver5.txt quiver6/combined-scores.quiver6.txt quiver7/combined-scores.quiver7.txt

for f in ../*united.txt; do b=`basename $f .txt`; echo $b; grep canu $f > $b.withextra.txt; grep falcon $f >> $b.withextra.txt; done
for f in *united.txt; do echo $f; b=`basename $f .txt`; mv $f $b.origrecipe.txt; done

cd ../canu_quiver_as_is
unitescores quiver1/combined-scores.quiver1.txt quiver2/combined-scores.quiver2.txt quiver3/combined-scores.quiver3.txt quiver4/combined-scores.quiver4.txt quiver5/combined-scores.quiver5.txt quiver6/combined-scores.quiver6.txt quiver7/combined-scores.quiver7.txt

for f in ../*united.txt; do b=`basename $f .txt`; echo $b; grep -E 'canu.default.pball.quiver|canu.minrl500.pball.quiver' $f > $b.withextra.txt; done
for f in *united.txt; do echo $f; b=`basename $f .txt`; mv $f $b.quivered_as_is.txt; done
for f in ../original_recipe/*origrecipe.txt; do b=`basename $f .txt`; echo $b; grep -E 'canu.default.pball.quiver|canu.minrl500.pball.quiver' $f > $b.origrecipe.txt; done

