for i in 4 11 13; do awk '{gsub(/,/,"\t"); print}' $1 | head -n1 | cut -f 1,$i
 grep -v ^# $1 | awk '{gsub(/,/,"\t"); print}'  | sort -k${i},${i}nr | cut -f 1,${i} | awk 'OFS="\t" {print $1,$2,NR}'; 
 echo
done
for i in 12; do 
 awk '{gsub(/,/,"\t"); print}' $1 | head -n1 | cut -f 1,$i 
 grep -v ^# $1 | awk '{gsub(/,/,"\t"); print}'  | sort -k${i},${i}n | cut -f 1,${i} | awk 'OFS="\t" {print $1,$2,NR}'
 echo
done
