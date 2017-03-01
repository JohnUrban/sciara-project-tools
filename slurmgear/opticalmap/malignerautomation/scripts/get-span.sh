for f in */merge/*Graph; do a=`echo $f | awk '{sub(/\/merge\//, "\t"); print $1}'`; s=`mergeBed -i $f | awk '{s+=$3-$2}END{print s}'`; echo -e $a"\t"$s; done
