for f in */reapr/output_directory/05.summary.report.txt; do 
s=`grep "Error free bases" $f | awk '{sub(/%/,""); print $4}'`; 
n=`echo $f | stringEdit - /reapr/output_directory/05.summary.report.txt`
echo -e $n "\t" $s; 
done
