for f in */ale/*.ALE.txt; do b=`basename $f .ALE.txt`; c=`head $f | grep ALE | awk '{print $3}'`; echo -e $b"\t"$c; done 
