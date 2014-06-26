#!/usr/bin/python

import foldx, re, shutil, random, os, math, sys, subprocess, glob
import numpy as np
from Bio import *
import Bio.PDB as PDB

def main():
  output = open('add_counts.txt', 'w')
  to_file = 'mutant\tabsolute_count\trelative_count\tbinding\tstab1\tstab2\n'
  output.write(to_file)
  output.close()

  all_tried_A = 0
  all_tried_C = 0
  all_tried_both = 0
  all_tried = {}

  #All of the tried mutants
  input_handle = open('all_mutants_tried.txt', 'r')
  for a_line in input_handle.readlines():
    split = a_line.split('\t')
    if split[0] == 'count':
      continue
    else:
      all_tried_both += 1

      if float((split[1])[1:-5]) <= 156:
        all_tried_A += 1
        if rmWS(split[1]) not in all_tried.keys():
          all_tried[rmWS(split[1])] = [all_tried_A, all_tried_both]
      else:
        all_tried_C += 1
        if rmWS(split[1]) not in all_tried.keys():
          all_tried[rmWS(split[1])] = [all_tried_C, all_tried_both]

  #All of the kept mutants
  input_handle = open('data.txt', 'r')
  output = open('add_counts.txt', 'a')
  for a_line in input_handle.readlines():
    split = a_line.split('\t')
    if split[0] == 'mutant':
      continue
    else:
      chain = ""
      if float((split[0])[1:-5]) <= 156:
        chain = '_A'
      else:
        chain = '_C'

      to_file = split[0] + chain +'\t' + str(all_tried[split[0]][1]) + '\t' + str(all_tried[split[0][0]) + '\t' + split[1] + '\t' + split[2] + '\t' + split[3] + '\n'
      output.write(to_file)
    
  output.close()

def rmWS(string):
    return("".join(string.split()))

#Run main program
if __name__ == '__main__':
   main()
