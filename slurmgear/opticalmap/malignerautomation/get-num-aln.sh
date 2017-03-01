for f in */merge/num*; do
  a=`echo $f | sed 's/\/merge\/num_alignments.txt//'`  
  n=`cat $f`
  echo -e $a"\t"$n
done
