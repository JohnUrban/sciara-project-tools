if [ $# -eq 0 ]; then SUF=comb; else SUF=$1; fi
for f in */sniffles_${SUF}*/*bedpe; do 
  c=`awk '{print $NF}' $f | grep -v pred | awkSum`; 
  n=`echo $f | awk '{sub(/\//,"\t"); print $1}'`
  echo -e $n"\t"$c; 
done 

#| sort -k2,2n
