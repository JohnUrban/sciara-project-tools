echo -e \#Rank"\t"Name"\t"NumComplete"\t"Total"\t"PctComplete
for f in */*/*/short_summary_* ; do n=`echo $f | awk '{sub(/\//,"\t"); print $1}'`; c=`grep "Complete BUSCOs" $f | awk '{print $1}'`; t=`grep "Total" $f | awk '{print $1}'`; echo -e $n"\t"$c"\t"$t; done | awk 'OFS="\t" {print $1,$2,$3,100.0*$2/$3}' | sort -k4,4nr | awk 'OFS="\t" {print NR,$0}' | sort -k2,2
