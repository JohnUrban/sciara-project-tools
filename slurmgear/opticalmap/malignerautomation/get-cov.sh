for f in */merge/*cov.txt; do
 a=`echo $f | sed 's/\/merge\/total_base_cov.txt//'`
 c=`cat $f`
 echo -e $a"\t"$c
done
