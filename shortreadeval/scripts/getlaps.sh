for file in */lap/*.lapscore; do 
 b=`cat $file`; 
 n=`basename $file .lapscore`
 echo -e $n"\t"$b | awk 'OFS="\t" {print $1,$2}'; 
done
