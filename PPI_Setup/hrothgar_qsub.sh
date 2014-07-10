#!/bin/csh

# Create Working Directory
set WDIR = $RESULTS/B_1K/B_1k_xxx
set FDIR = `pwd`

if ( -d $WDIR ) then
  rm -r $WDIR
endif

mkdir -p $WDIR

if ( ! -d $WDIR ) then
  echo $WDIR not created
  exit
endif

cd $WDIR

# Copy Data and Config Files
cp $FDIR/* .

echo '#\!/bin/bash' > hroth_run.q
echo '#$ -cwd' >> hroth_run.q
echo '#$ -V' >> hroth_run.q
echo '#$ -S /bin/bash' >> hroth_run.q
echo '#$ -N B_1k' >> hroth_run.q
echo '#$ -j y' >> hroth_run.q
echo '#$ -o main.log' >> hroth_run.q
echo '#$ -e main.err' >> hroth_run.q
echo '#$ -q normal' >> hroth_run.q
echo '#$ -pe fill 1' >> hroth_run.q
echo '#$ -P hrothgar\n' >> hroth_run.q

echo 'hostname\n' >> hroth_run.q

echo 'python main.py 2eke 1000 10 -23 -5.0 -9.7 final_data.txt\n' >> hroth_run.q

qsub hroth_run.q

