## NLINES should be 20 -- when 19, it is still running -- when >20, errors or cancellation present --- unfortunatel
## Should have 3 EDTs...

if [ $# -gt 0 ]; then 
  #wc -l $1/*/slu* | awk '$1!=20'
  grep -c EDT ${1}/*/sl* | awk '{sub(/:/,"\t"); print }' | awk '$2!=3 {OFS="\t"; gsub(/\//,"\t"); print $1"/"$2"/"$3,$2}' 
else
  #wc -l ipdsummary/*/slu* | awk '$1!=20'; 
  grep -c EDT ipdsummary/*/sl* | awk '{sub(/:/,"\t"); print }' | awk '$2!=3' | awk '$2!=3 {OFS="\t"; gsub(/\//,"\t"); print $1"/"$2"/"$3,$2}'
fi
