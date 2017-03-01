for f in */mreads/*err; do
 pct=`grep overall $f | awk '{sub(/%/,""); print $1}' `
 name=`echo $f | awk '{sub(/\//,"\t"); print $1}'`; 
 echo -e $name"\t"$pct
done 
