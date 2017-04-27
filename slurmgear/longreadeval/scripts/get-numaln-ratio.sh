#1=ont or pacbiot
for f in */snifflestats/$1.numaln; do
 a=`echo $f | sed 's/\/snifflestats\/.*numaln//g'`
 c=`cat $f`
 b=`basename $f .numaln`
 d=`cat $a/snifflestats/$1.numuniqaln`
 e=`echo $c/$d | bc -l`
 echo -e $a"\t"$c"\t"$d"\t"$e
done
