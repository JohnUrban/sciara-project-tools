for f in  */pilon/*changes; do name=`echo $f | awk '{sub(/\//,"\t"); print $1}'`; 
	pilonchanges.awk.sh $f | awk -v "name=$name" 'OFS="\t" {print name,$0}'
done > pilon2.changeinfo.txt
