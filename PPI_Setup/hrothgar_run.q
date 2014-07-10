#!/bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N B_1K_fold
#$ -j y
#$ -o main.log
#$ -e main.err
#$ -q normal
#$ -pe fill 1
#$ -P hrothgar
hostname

python main.py 2eke 1000 10 -23 -5.0 -9.7 final_data.txt
