#!/bin/bash

#define who and what here

MODEL=jules5.3

NAMELIST=`basename "$PWD"`

USER=ewanp82

#BSUB -o /work/scratch/${USER}/${MODEL}/${NAMELIST}/%J.o

#BSUB -e /work/scratch/${USER}/${MODEL}/${NAMELIST}/%J.e

#BSUB -q par-multi

#BSUB -n 1

#BSUB -W 1:00

# Get the original working directory

pwd1=${PWD}

# Get temp wd

pwd2=/work/scratch/${USER}/${MODEL}/${NAMELIST}

mkdir -p $pwd2

cp -dr * ${pwd2}/

cd $pwd2

#Define which JULES version to use

JULES=/home/users/${USER}/models/${MODEL}/build/bin/jules.exe

#Make some helpful output to screen and put the messages into files

start=`date`

echo "About to run 1 procs in \u2018pwd\u2018" >out.log

echo "Start Time: $start" >> out.log

echo /usr/local/bin/mpirun.lotus $JULES 2>> err.log >> out.log

module add parallel-netcdf/intel

/usr/local/bin/mpirun.lotus $JULES 2>> err.log >> out.log

end=`date`

echo "End Time: $end" >> out.log

echo "Storage MB: \u2018du -ms\u2018" >> out.log
