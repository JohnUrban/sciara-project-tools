function sub { echo $1 | sed "s/$2/$3$2/"; }

function sub0 ()
{ 
    new=`sub $1 $2 .quiver0x`
    mv $1 $new
}

to0x ()
{
    for d in $@; do
      if [ -d $d ]; then echo $d
          cd $d/ale
          sub0 *ale.err .ale.err
          sub0 *.ALE.txt .ALE.txt
          cd ../busco
          add0x run*
          cd run*
          add0x *
          cd ../../frc/
          sub0 *.frc_assemblyTable.csv .frc_assemblyTable.csv
          sub0 *.frcFeatures.gff .frcFeatures.gff
          cd ../lap
          sub0 *.lapscore .lapscore
          sub0 *.prob .prob
          cd ../mreads
          sub0 *.mapreads.err .mapreads.err
          cd ../../
       fi
     done
}

to0xauto ()
{
    to0x *quiver0x/
}

add0x () 
{ 
    for f in $@;
    do
      new=`echo $f | sed "s/\///g"`  
      mv $f $new.quiver0x;
    done
}


xto1x () 
{ 
    for f in $@;
    do
        new=`echo $f | sed 's/quiver/quiver1x/g'`;
        echo $new;
        mv $f $new;
    done
}



the4xALT () 
{ 
    for f in *ALT*;
    do
        new=`echo $f | sed 's/4xALT/quiver4x/g'`;
        echo $new;
        mv $f $new;
    done
}
the5xALT () 
{ 
    for f in *ALT*;
    do
        new=`echo $f | sed 's/5xALT/quiver5x/g'`;
        echo $new;
        mv $f $new;
    done
}
the6xALT () 
{ 
    for f in *ALT*;
    do
        new=`echo $f | sed 's/6xALT/quiver6x/g'`;
        echo $new;
        mv $f $new;
    done
}
the7xALT () 
{ 
    for f in *ALT*;
    do
        new=`echo $f | sed 's/7xALT/quiver7x/g'`;
        echo $new;
        mv $f $new;
    done
}
theALT () 
{ 
    for f in *ALT*;
    do
        new=`echo $f | sed 's/ALT//g'`;
        echo $new;
        mv $f $new;
    done
}
