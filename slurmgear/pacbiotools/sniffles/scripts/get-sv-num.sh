if [ $# -eq 0 ]; then SUF=comb; else SUF=$1; fi
for f in */sniffles_${SUF}*/*bedpe; do 
  c=`grep -c -v ^# $f `; 
  n=`echo $f | awk '{sub(/\//,"\t"); print $1}'`
  echo -e $n"\t"$c; 
done 

#| sort -k2,2n
