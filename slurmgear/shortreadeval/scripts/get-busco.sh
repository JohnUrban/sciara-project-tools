for f in */*/*/short_summary_* ; do 
 n=`echo $f | awk '{sub(/\//,"\t"); print $1}'`; 
 c=`grep "Complete BUSCOs" $f | awk '{print $1}'`; t=`grep "Total" $f | awk '{print $1}'`; 
 echo -e $n"\t"$c"\t"$t | awk 'OFS="\t" {print $1,100.0*$2/$3}' ; 
done 
