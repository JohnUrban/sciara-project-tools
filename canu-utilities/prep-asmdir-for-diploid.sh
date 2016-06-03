for file in g*/; do 
  echo $file; 
  rm $file/*fasta $file/*layout $file/unitigging.html; 
  rm -r $file/unitigging/*.tigStore/; 
  rm -r $file/unitigging/5-consensus/; 
  grep Guessed $file/unitigging/0-mercounts/*.ms22.estMerThresh.err; 
  echo; 
done


##make changes then
#for f in g*/; do echo $f; sh $f/unitigging/4-unitigger/unitigger.jobSubmit.sh; done
# then wait for jobs to finish and restart canu
