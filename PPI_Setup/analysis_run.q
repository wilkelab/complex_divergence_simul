#!/bin/bash
#$ -N WT_ancest
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
/share/apps/python-2.7.2/bin/python analyze_ppi.py final_data.txt 2eke.pdb

# Copy Results Back to Home Directory
mkdir -p $RDIR
cp * $RDIR

# Cleanup
rm -rf $WDIR
