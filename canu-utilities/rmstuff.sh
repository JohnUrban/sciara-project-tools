#!/bin/bash

## NEED TO COPY THIS INTO DIR ABOVE "sciara/"

diskUsageThisDirAndAllSubDir

echo "Top directory"
rm slurm*out &

echo "sciara"
cd sciara
rm *fastq *layout *Tig *Info &
rm -rf canu-logs/ &

echo "correction"
cd correction/
rm *err *kp *log *summary &

echo "correction/0-mercounts"
cd 0-mercounts/
rm *out *ignore *gram *info *gp *fasta &

echo "correction/1-overlapper"
cd ../1-overlapper/
rm *out *err *log &

echo "correction/2-correction"
cd ../2-correction/
rm *out *log *gp *summary *err *Correct &
rm -rf correction_inputs/ &
rm -rf correction_outputs/ &

echo "correction/sciara.gkpStore"
cd ../sciara.gkpStore/
rm blobs errorLog info* libraries load.dat *gp readNames.txt reads reads.txt &
cd ../
rm -rf sciara.ovlStore/ &

echo "trimming"
cd ../trimming/
rm *err *gkp *log *summary &
rm -rf sciara.ovlStore/ &

echo "trimming/0-mercounts"
cd 0-mercounts/
rm *out *err *fasta *gram *info *gp *ignore &

echo "trimming/1-overlapper"
cd ../1-overlapper/
rm *out *err *bat *job *opt &
rm -rf 001/ &


echo "trimming/3-overlapbasedtrimming"
cd ../3-overlapbasedtrimming/
rm * &

echo "trimming/sciara.gkpStore"
cd ../sciara.gkpStore/
rm blobs errorLog info* libraries load.dat *gp readNames.txt reads reads.txt &

echo "unitigging"
cd ../../unitigging/
rm *err *gkp *log *summary &
rm -rf sciara.ovlStore/ &

echo "unitigging/0-mercounts"
cd 0-mercounts/
rm *out *err *fasta *gram *info *gp &

echo "unitigging/1-overlapper/"
cd ../1-overlapper/
rm *out *err *bat *job *opt &
rm -rf 001/ &

echo "unitigging/3-overlapErrorAdjustment"
cd ../3-overlapErrorAdjustment/
rm *oea *out *err *.red *files &

echo "unitigging/4-unitigger"
cd ../4-unitigger/
rm *log *Info *map *ing *.ovl *err &

echo "unitigging/5-consensus"
echo "keeping conensus.out files"
cd ../5-consensus/
rm *files *err &

echo "unitigging/sciara.gkpStore"
cd ../sciara.gkpStore/
rm blobs errorLog info* libraries load.dat *gp readNames.txt reads* &
rm -rf partitions/ &

echo "unitigging/sciara.tigStore"
cd ../sciara.tigStore/
rm *dat *tig *filter &

cd ../../../
wait

diskUsageThisDirAndAllSubDir
