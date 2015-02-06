#!/bin/bash
#$ -N A
#$ -e /dev/null
#$ -o /dev/null
#$ -S /bin/bash

# Create Working Directory
WDIR=/state/partition1/$USER/$JOB_NAME-$JOB_ID
CDIR=/state/partition1/$USER
RDIR=/share/WilkeLab/work/agm854/Results/$JOB_NAME-$JOB_ID

FDIR=`pwd`

mkdir -p $WDIR

if [ ! -d $WDIR ]
then
  echo $WDIR not created
  exit
fi
cd $WDIR

# Copy Data and Config Files
cp $FDIR/* .

# Command to run
/share/apps/python-2.7.2/bin/python main.py 2eke random num_fixed no_selection 1000 10 10 -23.0 -5.0 -9.7 kept_mutants.txt all_mutants_tried.txt

# Copy Results Back to Home Directory
mkdir -p $RDIR
cp * $RDIR

# Cleanup
rm -rf $WDIR
