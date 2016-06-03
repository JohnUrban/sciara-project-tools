#!/bin/sh

jobid=$SLURM_ARRAY_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
  jobid=$1
fi
if [ x$jobid = x ]; then
  echo Error: I need SLURM_ARRAY_TASK_ID set, or a job index on the command line.
  exit 1
fi

if [ $jobid -gt 70 ]; then
  echo Error: Only 70 partitions, you asked for $jobid.
  exit 1
fi

jobid=`printf %04d $jobid`

syst=`uname -s`
arch=`uname -m`
name=`uname -n`

if [ "$arch" = "x86_64" ] ; then
  arch="amd64"
fi
if [ "$arch" = "Power Macintosh" ] ; then
  arch="ppc"
fi

bin="/gpfs_home/jurban/software/canu/canu/$syst-$arch/bin"

echo $syst-$arch-$name
echo $jobid
echo $SLURM_ARRAY_TASK_ID
