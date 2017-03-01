for f in */merge/score*; do
 c=`cat $f`
 a=`echo $f | sed 's/\/merge\/score.txt//'`
 echo -e $a"\t"$c
done
