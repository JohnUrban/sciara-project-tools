for f in */reapr/output_directory/per-base-mean-score.txt; do 
s=`cat $f`; 
n=`echo $f | stringEdit - /reapr/output_directory/per-base-mean-score.txt`; 
echo -e $n "\t" $s; 
done
