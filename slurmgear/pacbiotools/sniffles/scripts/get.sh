for f in */marginstats/$1; do
 a=`echo $f | sed 's/\/marginstats\/$1//'`
 c=`cat $f`
 echo -e $a"\t"$c
done
