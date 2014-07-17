#!/bin/bash
#$ -N B6_1k_9999B
#$ -e error_file
#$ -o out_file
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
/share/apps/python-2.7.2/bin/python main.py 2eke 1000 10 -23 -5.0 9999 final_data.txt 2> main.log 1> main.log

# Copy Results Back to Home Directory
mkdir -p $RDIR
cp * $RDIR

# Cleanup
rm -rf $WDIR
