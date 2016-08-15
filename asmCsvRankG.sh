for i in `awk '{gsub(/,/,"\t"); print}' $1 | cut -f 1 | grep -v ^# | paste -sd" "`; do 
bash asmCsvSummarizeG.sh $1 | grep ${i} | cut -f 3 | awkSum | awk -v "i=$i" '{print i"\t"$1}' 
done | sort -k2,2n
