for file in */frc/*gff; do
 if [ $# -eq 0 ] || [ $1 -eq 0 ]; then
#   a=`cat $file | wc -l`; b=`grep -v ^# $file | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; d=`echo $a/$g | bc -l`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$a"\t"$b"\t"$g"\t"$d"\t"$e
   b=`grep -v ^# $file | wc -l`; c=`basename $file .frcFeatures.gff`;
   #g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e
   echo -e $c"\t"$b
 elif [ $1 -eq 1 ]; then
  b=`grep -v ^# $file | awk '$3 == "COMPR_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 2 ]; then
  b=`grep -v ^# $file | awk '$3 == "HIGH_COV_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 3 ]; then
  b=`grep -v ^# $file | awk '$3 == "HIGH_NORM_COV_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 4 ]; then
  b=`grep -v ^# $file | awk '$3 == "HIGH_OUTIE_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 5 ]; then
  b=`grep -v ^# $file | awk '$3 == "HIGH_SINGLE_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 6 ]; then
  b=`grep -v ^# $file | awk '$3 == "HIGH_SPAN_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 7 ]; then
  b=`grep -v ^# $file | awk '$3 == "LOW_COV_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 8 ]; then
  b=`grep -v ^# $file | awk '$3 == "LOW_NORM_COV_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 elif [ $1 -eq 9 ]; then
  b=`grep -v ^# $file | awk '$3 == "STRECH_PE"' | wc -l`; c=`basename $file .frcFeatures.gff`; g=`cut -f 2 genomefiles/${c}.genome | awkSum`; e=`echo $b/$g | bc -l`; echo -e $c"\t"$b"\t"$g"\t"$e

 fi
done
